#include "diffbase.hpp"
#include <chrono>

#define DBL double

class Functor {
public:
    HOST_DEVICE DBL operator()(DBL x) const noexcept {
        return std::cos(2*x);
    }
};

HOST_DEVICE DBL func(DBL x) {
    return sin(10000*x);
}

HOST_DEVICE DBL func2(DBL x) {
    return x*x;
}

#ifdef __CUDACC__
template <typename F>
__global__ void kernel1(const F &f, double a, double b, uint64_t num) {
    DBL simpson = diff::integrate_simpson(f, a, b, num);
    printf("Simpson result: %lf\n", simpson);
}
template <typename F>
__global__ void kernel2(const F &f, double a, double b, double tol) {
    printf("Kernel 2 is running.\n");
    DBL trapq = diff::integrate_trapquad(f, a, b, tol);
    printf("Trapquad result: %lf\n", trapq);
}
#endif

int main(int argc, char **argv) {
    DBL a = -PI;
    DBL b =  PI;
    uint64_t num = 1'000'000;
    if (argc == 2)
        num = diff::misc::to_uint(*(argv + 1));
    Functor functor;
    std::cout << "With functor.\n" << std::endl;
    std::cout << "Integration range: [" << a << ',' << b << "]\n\n";
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    DBL trap = diff::integrate_trap(functor, a, b, num);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    std::chrono::duration trap_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    DBL tol = 1/256.0l;
    start = std::chrono::high_resolution_clock::now();
    DBL quad = diff::integrate_trapquad(functor, a, b, tol);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration quad_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    //num = 2;
    start = std::chrono::high_resolution_clock::now();
    DBL simpson = diff::integrate_simpson(functor, a, b, num);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration simp_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    uint64_t mdepth;
    DBL ptol = 1/65536.0;
    start = std::chrono::high_resolution_clock::now();
    DBL trapnr = diff::trapquad<double, Functor, 256>(functor, a, b, tol, ptol, &mdepth);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration nr_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Trapezium integration:\n\tNum. trapeziums = %" PRIu64 "\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\n"
           "Quad-trap integration:\n\tTolerance = %.20lf\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\n"
           "Non-recursive quad-trap:\n\tTolerance = %.20lf\n\tResult = %.20lf\n\tMax. depth = %" PRIu64
           "\n\tTime taken = %.20lf s\n\n",
           num, trap, trap_time.count()/((DBL) BILLION), tol, quad, quad_time.count()/((DBL) BILLION),
           tol, trapnr, mdepth, nr_time.count()/((double) BILLION));
    printf("Simpson's Rule:\n\tNum. parabolas = %" PRIu64 "\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\n",
        num, simpson, simp_time.count()/((DBL) BILLION));
#ifdef __CUDACC__
    dim3 blockSize{1, 1};
    dim3 gridSize{1, 1};
    kernel1<<<gridSize, blockSize>>>(functor, a, b, num);
    kernel2<<<gridSize, blockSize>>>(functor, a, b, tol);
    cudaError_t code;
    if ((code = cudaDeviceSynchronize()) != cudaSuccess) {
        fprintf(stderr, "Error: could not synchronise with device.\nReason: %s\n", cudaGetErrorString(code));
        return 1;
    }
    printf("After synchronisation.\n");
#endif
    short bin = 0b0110101001101001;
    char *s = diff::misc::binrepr(&bin, sizeof(short));
    std::cout << s << std::endl;
    delete [] s;
    return 0;
}
