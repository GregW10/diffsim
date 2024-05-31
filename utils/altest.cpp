#include "diffalloc.hpp"

int main(int argc, char **argv) {
    int num = -99999999;
    char buff[10];
    gtd::to_string(num, buff);
    std::cout << buff << std::endl;
    const char *s1 = "Hey there, bro!";
    const char *s2 = " What are you doing?\n";
    char both[256];
    gtd::strcpy_c(both, s1);
    gtd::strcat_c(both, s2);
    std::cout << both << std::endl;
    size_t count = 1024*1024;
    if (argc >= 2)
        count = diff::misc::to_uint(*(argv + 1));
    char *ptr = new char[count];
    int fd = open("shit.bin", O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR);
    gtd::write_all(fd, ptr, count);
    close(fd);
    delete [] ptr;
    diff::diff_alloc<long double> alloc{2000, 2000};
    std::cout << "Bytes written: " << alloc.to_dttr("test.dttr") << std::endl;
    return 0;
}
