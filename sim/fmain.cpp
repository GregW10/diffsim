#include "../utils/diffimg.hpp"

#define LAMBDA    0.000001l
#define YA       -0.000005l
#define YB        0.000005l
#define XA       -0.000005l
#define XB        0.000005l
#define DTTR_X    0.00002l
#define DTTR_Y    0.00002l
#define DTTR_DIST 0.00001l
#define DTTR_NX   10
#define DTTR_NY   10

#define LASER_POWER 0.005l // 5 mW laser
#define BEAM_DIAMETER 0.001l // 1 mm beam diameter

float gfunc(const float&) {
    return XA;
}

float hfunc(const float&) {
    return XB;
}

float beam_intensity(const float &laser_power, const float &beam_diam) {
    float r = beam_diam/2;
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
    float I0 = beam_intensity(LASER_POWER, BEAM_DIAMETER);
    printf("Incident light intensity: %.30f\n", I0);
    float abstol_y = 0.0625l/(1024.0l*1024.0l);
    float reltol_y = 0.0625l/(1024.0l*1024.0l);
    float ptol_y   = 1/(1024.0l*1024.0l*1024.0l);
    uint64_t mdepth_y    = diff::diffsim<float>::def_mdepth;
    float abstol_x = 0.0625l/(1024.0l*1024.0l);
    float reltol_x = 0.0625l/(1024.0l*1024.0l);
    float ptol_x   = 1/(1024.0l*1024.0l*1024.0l);
    uint64_t mdepth_x    = diff::diffsim<float>::def_mdepth;
    unsigned int threads = 1;
    unsigned int ptime   = 1;
    std::cout << "Starting simulation with the following parameters:"
                 "\n\tAperture lower x-limit = " << XA       << " m"
                 "\n\tAperture upper x-limit = " << XB       << " m"
                 "\n\tAperture lower y-limit = " << YA       << " m"
                 "\n\tAperture upper y-limit = " << YB       << " m"
                 "\n\tWavelength of light    = " << LAMBDA      << " m"
                 "\n\tDistance to detector   = " << DTTR_DIST        << " m"
                 "\n\tWidth of detector      = " << DTTR_X        << " m"
                 "\n\tLength of detector     = " << DTTR_Y        << " m"
                 "\n\tDetector x-resolution  = " << dttr_nx       <<
                 "\n\tDetector y-resolution  = " << dttr_ny       <<
                 "\n\tLight intensity        = " << I0       << " W/m^2"
                 "\n\tAbsolute y-tolerance   = " << abstol_y <<
                 "\n\tRelative y-tolerance   = " << reltol_y <<
                 "\n\tPeriodic y-tolerance   = " << ptol_y   <<
                 "\n\tAbsolute x-tolerance   = " << abstol_x <<
                 "\n\tRelative x-tolerance   = " << reltol_x <<
                 "\n\tPeriodic x-tolerance   = " << ptol_x   <<
                 "\n\tMax. y-recursion-depth = " << mdepth_y <<
                 "\n\tMax. x-recursion-depth = " << mdepth_x <<
                 "\n\tNumber of threads      = " << threads  <<
                 "\n\tProgress time delay    = " << ptime    << " s\n";
    diff::rectangle<float> ap{XA, XB, YA, YB};
    diff::diffimg<float> sim{LAMBDA, ap, DTTR_DIST, DTTR_X, DTTR_Y, I0, dttr_nx, dttr_ny};
    sim.diffract(abstol_y,//*1024.0l,
                 reltol_y,//*1024.0l,
                 ptol_y,
                 mdepth_y,
                 abstol_x,//*1024.0l,
                 reltol_x,//*1024.0l,
                 ptol_x,
                 mdepth_x,
                 threads,
                 ptime);
    // std::cout << "Max. depth reached in y: " << mdepth_y << ", in x: " << mdepth_x << std::endl;
    std::string path;
    sim.gen_bmp(path);
    std::cout << "BMP written to: " << path << std::endl;
    path.erase(path.end() - 3, path.end());
    path += "dffr";
    sim.to_dffr(path);
    std::cout << "DFFR written to: " << path << std::endl;
    gtd::complex<float> c{1, 1};
    // std::cout << gtd::abs(c) << std::endl;
    std::cout << c.mag() << std::endl;
    float s = 255.9;
    unsigned char x = s;
    std::cout << +x << std::endl;
    return 0;
}
