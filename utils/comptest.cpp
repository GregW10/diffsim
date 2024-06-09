#include "diffbase.hpp"
#include <iomanip>

gtd::complex<long double> func1(long double x, long double y) {
    return {std::cos(x), std::sin(y)};
}

gtd::complex<long double> func2(long double x, long double y) {
    gtd::complex<long double> num{x, y};
    return num*num;
}

long double gfunc(long double) {
    return -1;
}

long double hfunc(long double) {
    return 1;
}

int main() {
    gtd::complex<float> c{1, 1};
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.real(-1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.imag(-1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.real(1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c = c*5;
    std::cout << c << std::endl;
    c = 5*c;
    std::cout << c << std::endl;
    // c = {3, 2};
    gtd::complex<float> d{2, 3};
    std::cout << "(" << c << ")/(" << d << ") = " << c/d << std::endl;
    std::cout << "(" << c << ")*(" << c << ") = " << c*c << std::endl;
    long double ya = -1;
    long double yb = 1;
    gtd::complex<long double> res = diff::simpdblquad<long double, gtd::complex<long double>>(&func2, ya, yb, gfunc, hfunc);
    std::cout << std::setprecision(10);
    std::cout << res << std::endl;
    return 0;
}
