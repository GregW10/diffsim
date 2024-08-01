#include "diffsim.hpp"

static char buff[PATH_MAX + 8]{};

template <typename T> requires (std::is_floating_point_v<T>)
void print_info(const char *path) {
    std::cout << "File \"" << path << "\"\n";
    uint64_t slen = gtd::strlen_c(path) + 7;
    gtd::write_all(STDOUT_FILENO, buff, slen);
    std::cout << diff::diffsim<T>{path, true, false} << std::endl;
}

int main(int argc, char **argv) {
    uint64_t counter = PATH_MAX + 7;
    char *ptr = buff;
    while (counter --> 0)
        *ptr++ = '-';
    while (--argc > 0) {
        try {
            try {
                print_info<long double>(*argv);
            } catch (const diff::invalid_dffr_format&) {
                try {
                    print_info<double>(*argv);
                } catch (const diff::invalid_dffr_format&) {
                    try {
                        print_info<float>(*argv);
                    } catch (const diff::invalid_dffr_format &_e) {
                        std::cerr << "Error: file \"" << *argv << "\" could not be converted due to invalid file "
                                                                  "format.\nReason:\n" << _e.what() << '\n';
                    }
                }
            }
        } catch (const std::exception &_e) {
            std::cerr << "Error: file \"" << *argv << "\" could not be converted.\nReason:\n" << _e.what() << '\n';
        }
    }
    return 0;
}
