#include "diffbase.hpp"
#include <chrono>

long double func(const long double &x) {
    return sinl(10000*x);
}

long double func2(const long double &x) {
    return x*x;
}

int main(int argc, char **argv) {
    long double a = 0;
    long double b = 3;
    uint64_t num = 1'000'000;
    if (argc == 2)
        num = diff::misc::to_uint(*(argv + 1));
    std::cout << "Integration range: [" << a << ',' << b << "]\n\n";
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    long double trap = diff::integrate_trap(&func, a, b, num);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    std::chrono::duration trap_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    long double tol = 0.0001l;
    start = std::chrono::high_resolution_clock::now();
    long double quad = diff::integrate_trapquad(&func, a, b, tol);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration quad_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    num = 2;
    start = std::chrono::high_resolution_clock::now();
    long double simpson = diff::integrate_simpson(&func2, a, b, num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration simp_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Trapezium integration:\n\tNum. trapeziums = %" PRIu64 "\n\tResult = %Lf\n\tTime taken = %Lf s\n\nQuad-trap "
           "integration:\n\tTolerance = %Lf\n\tResult = %Lf\n\tTime taken = %Lf s\n\n",
           num, trap, trap_time.count()/((long double) BILLION), tol, quad, quad_time.count()/((long double) BILLION));
    printf("Simpson's Rule:\n\tNum. parabolas = %" PRIu64 "\n\tResult = %Lf\n\tTime taken = %Lf s\n\n",
        num, simpson, simp_time.count()/((long double) BILLION));
    return 0;
}
