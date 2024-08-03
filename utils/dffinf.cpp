#include "diffsim.hpp"

static char buff[PATH_MAX + 8]{};

template <typename T> requires (std::is_floating_point_v<T>)
void print_info(const char *path) {
    if constexpr (std::same_as<T, long double>) {
        std::cout << "File \"" << path << "\"\n";
        uint64_t slen = gtd::strlen_c(path) + 7;
        gtd::write_all(STDOUT_FILENO, buff, slen);
        putchar('\n');
    }
    std::cout << diff::diffsim<T>{path, true, false} << std::endl;
}

int main(int argc, char **argv) {
    uint64_t counter = PATH_MAX + 7;
    char *ptr = buff;
    while (counter --> 0)
        *ptr++ = '-';
    ++argv;
    uint64_t sread = 0;
    uint64_t niread = 0;
    uint64_t oiread = 0;
    while (--argc > 0) {
        try {
            try {
                print_info<long double>(*argv);
                ++sread;
            } catch (const diff::invalid_dffr_format&) {
                try {
                    print_info<double>(*argv);
                    ++sread;
                } catch (const diff::invalid_dffr_format&) {
                    try {
                        print_info<float>(*argv);
                        ++sread;
                    } catch (const diff::invalid_dffr_format &_e) {
                        std::cerr << "Error: file \"" << *argv << "\" could not be converted due to invalid file "
                                                                  "format.\nReason:\n" << _e.what() << '\n';
                        ++niread;
                    }
                }
            }
        } catch (const std::exception &_e) {
            std::cerr << "Error: file \"" << *argv << "\" could not be converted.\nReason:\n" << _e.what() << '\n';
            ++oiread;
        }
        ++argv;
    }
    printf("\033[1m\033[38;5;10mFiles successfully read: \033[1m\033[38;5;11m%" PRIu64 "\n"
           "\033[38;5;9mFiles failed to be read due to\n\t"
           "\033[38;5;14mInvalid .dffr format: \033[1m\033[38;5;11m%" PRIu64 "\n\t"
           "\033[38;5;214mIO exception:         \033[1m\033[38;5;11m%" PRIu64 "\n\033[0m", sread, niread, oiread);
    return niread ? (oiread ? 255 : 127) : (oiread ? 1 : 0);
}
