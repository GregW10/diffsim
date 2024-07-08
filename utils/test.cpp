#include "diffbase.hpp"
#include "../glib/misc/gregstack.hpp"
#include <chrono>
#include <complex>

#define DBL double

class Functor {
    mutable uint64_t times_called = 0;
public:
    HOST_DEVICE DBL operator()(DBL x) const noexcept {
        ++times_called;
        DBL val = x;
        return std::cos(x);
    }
    uint64_t reset() const noexcept {
        uint64_t val = times_called;
        times_called = 0;
        return val;
    }
    uint64_t ncalls() const noexcept {
        return times_called;
    }
};

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
    DBL a = 0.000;
    DBL b = 16*PI;
    uint64_t num = 10'000'000;
    uint64_t mdepth = 32;
    // char *endptr;
    if (argc == 2)
        mdepth = diff::misc::to_uint(*(argv + 1));
    uint64_t mdepth_s = mdepth;
    // if (endptr == *(argv + 1))
    //     return 1;
    Functor functor;
    DBL abstol = (1/65536.0l);///16384.0l;///65536.0l;
    DBL reltol = 0.0000001l;
    DBL ptol = 1/(1024.0l*1024.0l);
    std::cout << "Integration range: [" << a << ',' << b << "]\n\n";
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    DBL trap = diff::integrate_trap(functor, a, b, num);
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    uint64_t fc_trap = functor.reset();
    std::chrono::duration trap_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    start = std::chrono::high_resolution_clock::now();
    DBL quad = 1;//diff::integrate_trapquad(functor, a, b, tol, ptol);
    end = std::chrono::high_resolution_clock::now();
    uint64_t fc_trapr = functor.reset();
    std::chrono::duration quad_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    //num = 2;
    start = std::chrono::high_resolution_clock::now();
    DBL simpson = diff::integrate_simpson(functor, a, b, num);
    end = std::chrono::high_resolution_clock::now();
    uint64_t fc_simps = functor.reset();
    std::chrono::duration simp_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // functor.reset();
    /*std::ofstream out{"intervals.bin", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc};*/
    start = std::chrono::high_resolution_clock::now();
    DBL trapnr = 5;//diff::trapquad(functor, a, b, abstol, ptol, &mdepth/*, &out*/);
    end = std::chrono::high_resolution_clock::now();
    /*out.close();*/
    uint64_t fc_trapi = functor.reset();
    std::chrono::duration nr_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    start = std::chrono::high_resolution_clock::now();
    DBL simpr = diff::simpquad<DBL, DBL, decltype(functor), true>(functor, a, b, abstol, reltol, ptol, &mdepth_s/*, &out*/);
    end = std::chrono::high_resolution_clock::now();
    /*out.close();*/
    uint64_t fc_simpr = functor.reset();
    std::chrono::duration simpr_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Trapezium integration:\n\tNum. trapeziums = %" PRIu64 "\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\t"
                                                                  "Num. func. calls = %" PRIu64 "\n\n"
           "Quad-trap integration:\n\tTolerance = %.20lf\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\t"
           "Num. func. calls = %" PRIu64 "\n\n"
           "Non-recursive quad-trap:\n\tTolerance = %.20lf\n\tResult = %.20lf\n\tMax. depth = %" PRIu64
           "\n\tTime taken = %.20lf s\n\tNum. func. calls = %" PRIu64 "\n\n",
           num, trap, trap_time.count()/((DBL) BILLION), fc_trap, abstol, quad, quad_time.count()/((DBL) BILLION),
           fc_trapr, abstol, trapnr, mdepth, nr_time.count()/((double) BILLION), fc_trapi);
    printf("Simpson's Rule:\n\tNum. parabolas = %" PRIu64 "\n\tResult = %.20lf\n\tTime taken = %.20lf s\n\t"
                                                          "Num. func. calls = %" PRIu64 "\n\n",
        num, simpson, simp_time.count()/((DBL) BILLION), fc_simps);
    printf("Simpson's AQ Rule:\n\tAbs. Tolerance = %.20lf\n\tRel. Tolerance = %.20lf\n\tResult = %.20lf\n\t"
           "Max. depth = %" PRIu64 "\n\tTime taken = %.20lf s\n\tNum. func. calls = %" PRIu64 "\n\n",
        abstol, reltol, simpr, mdepth_s, simpr_time.count()/((DBL) BILLION), fc_simpr);
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
