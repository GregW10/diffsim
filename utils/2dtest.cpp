#include "diffbase.hpp"
#include <chrono>

#define DBL double

#define DFORMAT "lf"

class Functor2D {
    mutable uint64_t times_called = 0;
public:
    HOST_DEVICE DBL operator()(DBL x, DBL y) const noexcept {
        ++times_called;
        return std::sin(x)*std::cos(y);
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

DBL gfunc(DBL y) {
    y -= 1;
    //printf("gfunc: %lf\n", 1 - std::sqrt(1 - y*y));
    return 1 - std::sqrt(1 - y*y);
}

DBL hfunc(DBL y) {
    y -= 1;
    //printf("hfunc: %lf\n", 1 + std::sqrt(1 - y*y));
    return 1 + std::sqrt(1 - y*y);
}

int main(int argc, char **argv) {
    uint64_t mdepth_x = 36;
    uint64_t mdepth_y = 36;
    if (argc == 2)
        mdepth_x = mdepth_y = diff::misc::to_uint(*(argv + 1));
    else if (argc == 3) {
        mdepth_x = diff::misc::to_uint(*(argv + 1));
        mdepth_y = diff::misc::to_uint(*(argv + 2));
    }
    else if (argc > 3)
        return 1;
    Functor2D f;
    DBL ya = 0;
    DBL yb = 2;
    DBL abstol_y = 0.0625l/102400000.0l;
    DBL reltol_y = 0.0625l/102400000.0l;
    DBL ptol_y = 0.00000000000000000001;
    DBL abstol_x = 0.0625l/10240000000000000.0l;
    DBL reltol_x = 0.0625l/10240000000000000.0l;
    DBL ptol_x = 0.000000000000000000000000001;
    std::chrono::time_point<std::chrono::high_resolution_clock> start =
            std::chrono::high_resolution_clock::now();
    DBL res = diff::simpdblquad<DBL, DBL, Functor2D, decltype(gfunc), true>
            (f, ya, yb, gfunc, hfunc, abstol_y, reltol_y, ptol_y, &mdepth_y, abstol_x, reltol_x, ptol_x, &mdepth_x);
    std::chrono::time_point<std::chrono::high_resolution_clock> end =
            std::chrono::high_resolution_clock::now();
    std::chrono::duration tot = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    printf("Result: %.20" DFORMAT "\n", res);
    printf("Num func. calls: %" PRIu64 "\n", f.ncalls());
    printf("Max. depth in y: %" PRIu64 "\n", mdepth_y);
    printf("Max. depth in x: %" PRIu64 "\n", mdepth_x);
    std::cout << "Time taken: " << tot.count()/((long double) BILLION) << " s" << std::endl;
    return 0;
}
