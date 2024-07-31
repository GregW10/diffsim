#include "diffsim.hpp"
#include <set>
#include <random>

#define DEF_RANDLEN 10

uint64_t randlen = DEF_RANDLEN;
bool only_dffr = false; // determines whether to only check files ending in ".dffr" (in case they don't have that ext.)

std::string *rand_fname(std::string *out) {
    if (!out)
        return nullptr;
    static std::set<std::string> taken_names;
    std::string rand_name;
    rand_name.reserve(randlen + 5);
    // generate random mid-part
    std::mt19937_64 mersenne{std::random_device{}()};
    std::uniform_int_distribution<char> lletters{97, 122};
    std::uniform_int_distribution<char> uletters{65, 90};
    std::uniform_int_distribution<char> digits{48, 57};
    std::uniform_int_distribution<char> three{1, 3};
    uint64_t counter;
    unsigned char _t;
    auto end_it = taken_names.end();
    do {
        counter = randlen;
        while (counter --> 0) {
            if ((_t = three(mersenne)) == 1)
                rand_name += lletters(mersenne);
            else if (_t == 2)
                rand_name += uletters(mersenne);
            else
                rand_name += digits(mersenne);
        }
    } while(taken_names.find(rand_name) != end_it);
    taken_names.insert(rand_name);
    out->operator+=(rand_name);
    out->operator+=(".dffr");
    return out;
}

const char *get_dirname(int ac, char **av, int *indices) {
    const char *dname{};
    while (ac --> 0) {
        if (**av == '-') {
            if (*(*av + 1) == 'o') {
                if (!ac) {
                    fprintf(stderr, "Error: \"-o\" flag passed without subsequent argument.\n");
                    exit(1);
                }
                // *endptr = av;
                *indices = ac;
                dname = *++av;
                --ac;
            }
            else if (*(*av + 1) == 'l') {
                if (!ac) {
                    fprintf(stderr, "Error: \"-l\" flag passed without subsequent argument.\n");
                    exit(1);
                }
                uint64_t num_digs = (uint64_t) (log10l(NAME_MAX - 5) + 1.0l);
                if (gtd::strlen_c(*++av) > num_digs) {
                    fprintf(stderr, "Error: two many digits in random-name-length passed.\n");
                    exit(1);
                }
                if (!gtd::to_uint(*av, &randlen)) {
                    fprintf(stderr, "Error: could not convert random-name-length argument \"%s\" to unsigned integral "
                                    "form.\n", *av);
                    exit(1);
                }
                if (randlen > NAME_MAX) {
                    fprintf(stderr, "Error: random-name-length of %" PRIu64 " exceeds system limit of %" PRIu64 ".\n",
                        randlen, (uint64_t) NAME_MAX);
                    exit(1);
                }
                if (!randlen) {
                    fprintf(stderr, "Error: random-name-length cannot be zero.\n");
                    exit(1);
                }
                *(indices + 1) = ac;
                --ac;
            }
            else if (*(*av + 1) == 'd') {
                only_dffr = true;
                *(indices + 2) = ac;
            }
        }
        ++av;
    }
    if (!dname) {
        // *endptr = nullptr;
        static char cbuff[26 + 7 + 8]{};
        dname = cbuff;
        gtd::strcpy_c(cbuff, gtd::get_time());
        gtd::strcat_c(cbuff, "_fconv_");
        gtd::strcat_c(cbuff, "XXXXXXXX"); // results in 218'340'105'584'896 different possible filenames
        if (!mkdtemp(cbuff)) {
            fprintf(stderr, "Error: could not create output directory \"%s\".\n", cbuff);
            exit(1);
        }
    } else {
        struct stat buff{};
        if (stat(dname, &buff) == -1) {
            fprintf(stderr, "Error: could not obtain file information for \"%s\".\n", dname);
            exit(1);
        }
        if (!S_ISDIR(buff.st_mode)) {
            fprintf(stderr, "Error: file \"%s\" is NOT a directory.\n", dname);
            exit(1);
        }
    }
    return dname;
}

template <typename T> requires (std::is_floating_point_v<T>)
bool convert_directory(const char*, const char*, uint64_t*, bool = false);

template <typename T> requires (std::is_floating_point_v<T>)
void convert_file(const char *from, const char *ddir, uint64_t *convf) {
    struct stat buff{};
    if (stat(from, &buff) == -1) {
        fprintf(stderr, "Error: could not obtain \"%s\" file information.\n", from);
        return;
    }
    if (S_ISDIR(buff.st_mode)) {
        convert_directory<T>(from, ddir, convf, false);
        return;
    }
    if (!S_ISREG(buff.st_mode)) {
        fprintf(stderr, "Error: file \"%s\" is not a regular file.\n", from);
        return;
    }
    std::string opath = ddir;
    if (!opath.ends_with('/'))
        opath += '/';
    // opath += from;
    rand_fname(&opath);
    std::cout << "Attempting to create: " << opath << std::endl;
    std::cout << "Attempting to read: " << from << std::endl;
    try {
        try {
            diff::diffsim<long double> insim{from, false};
            diff::diffsim<T>{insim}.to_dffr(opath);
        } catch (const diff::invalid_dffr_format&) {
            try {
                diff::diffsim<double> insim{from, false};
                diff::diffsim<T>{insim}.to_dffr(opath);
            } catch (const diff::invalid_dffr_format&) {
                try {
                    diff::diffsim<float> insim{from, false};
                    diff::diffsim<T>{insim}.to_dffr(opath);
                } catch (const diff::invalid_dffr_format &_e) {
                    std::cerr << "Error: file \"" << from << "\" could not be converted. Reason:\n"<< _e.what() << '\n';
                    return;
                }
            }
        }
    } catch (const std::exception &_e) {
        std::cerr << "Error: file \"" << from << "\" could not be converted. Reason:\n" << _e.what() << '\n';
        return;
    }
    ++*convf;
}

template <typename T> requires (std::is_floating_point_v<T>)
bool convert_directory(const char *dirname, const char *ddir, uint64_t *convf, bool have_ddir) {
    struct dirent *entry;
    DIR *dir = opendir(dirname);
    if (!dir) {
        fprintf(stderr, "Error: could not open directory \"%s\".\n", dirname);
        return false;
    }
    // GETCWD AND THEN RETURN TO THAT AT THE END OF THE FUNCTION!!!!
    if (chdir(dirname) == -1) {
        closedir(dir);
        fprintf(stderr, "Error: could not change working directory to \"%s\".\n", dirname);
        return false;
    }
    bool have_dot = false;
    bool have_dot_dot = false;
    // bool have_ddir = false || have_ddir_parent;
    while ((entry = readdir(dir))) {
        if (!have_dot && gtd::str_eq(entry->d_name, ".")) {
            have_dot = true;
            continue;
        }
        if (!have_dot_dot && gtd::str_eq(entry->d_name, "..")) {
            have_dot_dot = true;
            continue;
        }
        if (entry->d_type == DT_DIR) {
            if (!have_ddir) {
                char buff[PATH_MAX]{};
                if (realpath(entry->d_name, buff) == nullptr) {
                    fprintf(stderr, "Error: could not resolve pathname for directory \"%s\".\n", entry->d_name);
                    exit(1);
                }
                if (gtd::str_eq(buff, ddir)) {
                    have_ddir = true;
                    continue;
                }
            }
            // if (stat(entry->d_name, &buff) == -1) {
            //     fprintf(stderr, "Error: could not obtain file information for \"%s\".\n", entry->d_name);
            //     continue;
            // }
            if (convert_directory<T>(entry->d_name, ddir, convf, have_ddir)) {
                if (chdir("..") == -1) {
                    fprintf(stderr, "Error: could not move up to parent directory.\n");
                    exit(1);
                }
            }
        } else
            convert_file<T>(entry->d_name, ddir, convf);
    }
    closedir(dir);
    return true;
}

template <typename T> requires (std::is_floating_point_v<T>)
void convert_files(int argc, char **argv) {
    argc -= 2;
    argv += 2;
    // char **endptr{};
    int indices[3] = {argc, argc, argc};
    const char *ddir = get_dirname(argc, argv, indices);
    char ddir_full[PATH_MAX]{};
    if (realpath(ddir, ddir_full) == nullptr) {
        fprintf(stderr, "Error: could not obtain fully resolved pathname for target directory \"%s\".\n", ddir);
        exit(1);
    }
    uint64_t convf = 0;
    while (argc --> 0) {
        // std::cout << "*argv: " << *argv << std::endl;
        if (argc == indices[0] || argc == indices[1]) {
            argv += 2;
            --argc;
            continue;
        }
        if (argc == indices[2]) {
            ++argv;
            continue;
        }
        convert_file<T>(*argv++, ddir_full, &convf);
    } /*
    if (endptr) {
        while (argc --> 0) {
            if (endptr == argv) {
                argv += 2;
                --argc;
                continue;
            }
            convert_file<T>(*argv++, ddir_full, &convf);
        }
    }
    else
        while (argc --> 0)
            convert_file<T>(*argv++, ddir_full, &convf); */
    printf("Total number of files converted: %" PRIu64 "\n", convf);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Error: insufficient arguments.\n"
                        "Usage: ./fconv <to> [-o <DIR>] <dffr_file1> [<more_dffr_files>...]\n");
        return 1;
    }
    const char *to = *(argv + 1);
    if (gtd::str_eq(to, "f"))
        convert_files<float>(argc, argv);
    else if (gtd::str_eq(to, "lf"))
        convert_files<double>(argc, argv);
    else if (gtd::str_eq(to, "Lf"))
        convert_files<long double>(argc, argv);
    else {
        fprintf(stderr, "Error: unrecognised floating-point type \"%s\".\n"
                        "Options are:\n\tfloat: \"f\"\n\tdouble: \"lf\"\n\tlong double: \"Lf\"\n", to);
        return 1;
    }
    return 0;
}
