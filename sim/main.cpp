#include "../utils/diffimg.hpp"

#define LAMBDA    0.000001l
#define YA       -0.000005l
#define YB        0.000005l
#define XA       -0.000005l
#define XB        0.000005l
#define DTTR_X    0.00002l
#define DTTR_Y    0.00002l
#define DTTR_DIST 0.00001l
#define DTTR_NX 10
#define DTTR_NY 10

#define LASER_POWER 0.005l // 5 mW laser
#define BEAM_DIAMETER 0.001l // 1 mm beam diameter

long double gfunc(const long double&) {
    return XA;
}

long double hfunc(const long double&) {
    return XB;
}

long double beam_intensity(const long double &laser_power, const long double &beam_diam) {
    long double r = beam_diam/2;
    return laser_power/(PI*r*r);
}

int main(int argc, char **argv) {
    uint64_t dttr_nx = DTTR_NX;
    uint64_t dttr_ny = DTTR_NY;
    if (argc == 3) {
        dttr_nx = diff::misc::to_uint(*(argv + 1));
        dttr_ny = diff::misc::to_uint(*(argv + 2));
    } else if (argc != 1)
        return 1;
    std::cout << "DTTR_NX: " << dttr_nx << ", DTTR_NY: " << dttr_ny << std::endl;
    long double I0 = beam_intensity(LASER_POWER, BEAM_DIAMETER);
    diff::aperture<long double> ap{YA, YB, &gfunc, &hfunc};
    diff::diffimg<long double> sim{LAMBDA, ap, DTTR_DIST, DTTR_X, DTTR_Y, I0, dttr_nx, dttr_ny};
    sim.diffract();
    sim.gen_bmp();
    gtd::complex<long double> c{1, 1};
    // std::cout << gtd::abs(c) << std::endl;
    std::cout << c.mag() << std::endl;
    long double s = 255.9;
    unsigned char x = s;
    std::cout << +x << std::endl;
    return 0;
}
