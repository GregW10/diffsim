#include "../utils/diffbase.hpp"

double xa = 1;
double xb = 3;

gtd::complex<double> func(double x) {
    return {std::log(x), std::cbrt(x)*std::cos(x)};
}

gtd::complex<double> func2(double x, double y) {
    return {y*std::sqrt(y) + std::log(x), std::sqrt(x*y)*std::exp(x*x*y*y)};
}

double gfunc(const double&) {
    return xa;
}

double hfunc(const double&) {
    return xb;
}

int main() {
    double a = 1;
    double b = 3;
    uint64_t num = 20'000'000;
    gtd::complex<double> ret = diff::simpson<double, gtd::complex<double>, decltype(&func)>(&func, a, b, num);
    std::cout << "Result: " << ret << std::endl;
    double ya = 2;
    double yb = 4;
    uint64_t num_y = 20'000;
    uint64_t num_x = 20'000;
    gtd::complex<double> ret2 =
        diff::dbl_simpson<double, gtd::complex<double>, decltype(&func2), decltype(gfunc), decltype(hfunc)>
            (func2, ya, yb, gfunc, hfunc, num_y, num_x);
    std::cout << "Result: " << ret2 << std::endl;
    return 0;
}
