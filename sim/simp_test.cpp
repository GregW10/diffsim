#include "../utils/diffbase.hpp"
#include <chrono>

double xa = 1;
double xb = 3;

gtd::complex<double> func(double x) {
    return {std::cos(x), std::sin(x)};
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
    std::cout.precision(20);
    double a = 0;
    double b = 128*PI;
    uint64_t num = 2'000'000;
    gtd::complex<double> ret = diff::simpson<double, gtd::complex<double>, decltype(&func)>(&func, a, b, num);
    std::cout << "Result: " << ret << std::endl;
    double abstol_x = 1e-12;
    double reltol_x = 1e-12;
    double ptol_x = 1e-9;
    uint64_t mdepth_x = 32;
    std::chrono::time_point<std::chrono::high_resolution_clock> start =
        std::chrono::high_resolution_clock::now();
    gtd::complex<double> r = diff::simpqr<double, gtd::complex<double>, decltype(func)>
        (func, a, b, abstol_x, reltol_x, ptol_x, &mdepth_x);
    std::chrono::time_point<std::chrono::high_resolution_clock> end =
        std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << "Recursive result: " << r << "\nTime taken: " <<
        ((long double) elapsed.count())/BILLION << " s\nMax. depth: " << mdepth_x << std::endl;
    double ya = 2;
    double yb = 4;
    uint64_t num_y = 2'000;
    uint64_t num_x = 2'000;
    gtd::complex<double> ret2 =
        diff::dbl_simpson<double, gtd::complex<double>, decltype(&func2), decltype(gfunc), decltype(hfunc)>
            (func2, ya, yb, gfunc, hfunc, num_y, num_x);
    std::cout << "Result: " << ret2 << std::endl;
    double abstol_y = 1e-12;
    double reltol_y = 1e-12;
    double ptol_y = 1e-9;
    abstol_x = abstol_y;
    reltol_x = reltol_y;
    ptol_x = ptol_y;
    uint64_t mdepth_y = 32;
    mdepth_x = mdepth_y;
    start = std::chrono::high_resolution_clock::now();
    gtd::complex<double> r2 =
        diff::simpdqr<double, gtd::complex<double>, decltype(&func2), decltype(gfunc), decltype(hfunc)>
            (func2, ya, yb, gfunc, hfunc, abstol_y, reltol_y, ptol_y, &mdepth_y, abstol_x, reltol_x, ptol_x, &mdepth_x);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << "Recursive result: " << r2 << "\nMax. y-depth: " << mdepth_y << "\nMax. x-depth: " << mdepth_x <<
        "\nTime taken: " << ((long double)elapsed.count())/BILLION << " s" << std::endl;
    return 0;
}
