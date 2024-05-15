#include "diffbase.hpp"

long double func(const long double &x) {
    return x*x*x;
}

int main(int argc, char **argv) {
    long double a = -3;
    long double b = 3;
    uint64_t num = 1'000'000'000;
    if (argc == 2)
        num = diff::misc::to_uint(*(argv + 1));
    std::cout << "The integral of \"func\" from x = " << a << " to x = " << b << " with " << num << " trapeziums is: " <<
    diff::integrate_trap(&func, a, b, num) << std::endl;
    return 0;
}
