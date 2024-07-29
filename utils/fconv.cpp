#include "diffsim.hpp"

const char *get_dirname(int ac, char **av) {
    const char *dname{};
    while (ac --> 0) {
        if (**av == '-' && *(*av + 1) == 'o') {
            if (!ac) {
                fprintf(stderr, "Error: \"-o\" flag passed without subsequent argument.\n");
                exit(1);
            }
            dname = *++av;
            --ac;
        }
        ++av;
    }
    if (!dname) {
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
void convert_directory(const char*, uint64_t*);

template <typename T> requires (std::is_floating_point_v<T>)
void convert_file(const char *from, const char *ddir, uint64_t *convf) {
    struct stat buff{};
    if (stat(from, &buff) == -1) {
        fprintf(stderr, "Error: could not obtain \"%s\" file information.\n", from);
    }
    if (S_ISDIR(buff.st_mode)) {
        convert_directory<T>(from, convf);
    }
    if (!S_ISREG(buff.st_mode)) {
        fprintf(stderr, "Error: file \"%s\" is not a regular file.\n", from);
    }

    ++*convf;
}

template <typename T> requires (std::is_floating_point_v<T>)
void convert_directory(const char *dir, uint64_t *convf) {

}

template <typename T> requires (std::is_floating_point_v<T>)
void convert_files(int argc, char **argv) {
    argc -= 2;
    argv += 2;
    const char *ddir = get_dirname(argc, argv);
    uint64_t convf = 0;
    while (argc --> 0) {
        convert_file<T>(*argv++, ddir, &convf);
    }
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
