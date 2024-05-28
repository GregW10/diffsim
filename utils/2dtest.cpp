#include "diffbase.hpp"

#define DBL double

class Functor2D {
    mutable uint64_t times_called = 0;
public:
    HOST_DEVICE DBL operator()(DBL x, DBL y) const noexcept {
        ++times_called;
        x -= 1;
        y -= 1;
        return x*x + y*y*y*y;
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
    uint64_t mdepth_x = NO_MAX_DEPTH;
    uint64_t mdepth_y = NO_MAX_DEPTH;
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
    DBL abstol_y = 0.0625l/10240000.0l;
    DBL reltol_y = 0.0625l/10240000.0l;
    DBL ptol_y = -0.000000000000001;
    DBL abstol_x = 0.0625l/10240000.0l;
    DBL reltol_x = 0.0625l/10240000.0l;
    DBL ptol_x = -0.000000000000001;
    DBL res = diff::simpdblquad(f, ya, yb, gfunc, hfunc, abstol_y, reltol_y, ptol_y, &mdepth_y,
                                abstol_x, reltol_x, ptol_x, &mdepth_x);
    printf("Result: %.20lf\n", res);
    printf("Num func. calls: %" PRIu64 "\n", f.ncalls());
    printf("Max. depth in y: %" PRIu64 "\n", mdepth_y);
    printf("Max. depth in x: %" PRIu64 "\n", mdepth_x);
    return 0;
}
