#include "../utils/diffimg.hpp"

#define LAMBDA    0.000001l
#define YA       -0.000005l
#define YB        0.000005l
#define XA       -0.000005l
#define XB        0.000005l
#define DTTR_X    0.00002l
#define DTTR_Y    0.00002l
#define DTTR_DIST 0.00001l
#define DTTR_NX 100
#define DTTR_NY 100

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

int main() {
    std::cout << "DTTR_NX: " << DTTR_NX << ", DTTR_NY: " << DTTR_NY << std::endl;
    long double I0 = beam_intensity(LASER_POWER, BEAM_DIAMETER);
    diff::aperture<long double> ap{YA, YB, &gfunc, &hfunc};
    diff::diffimg<long double> sim{LAMBDA, ap, DTTR_DIST, DTTR_X, DTTR_Y, I0, DTTR_NX, DTTR_NY};
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
