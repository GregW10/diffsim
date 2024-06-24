#include "../utils/diffimg.hpp"
#include "../glib/misc/gregparse.hpp"

#define DEF_NX   500
#define DEF_NY   500
#define DEF_LAM  0.000001l // blue-ish purple light
#define DEF_XA  -0.000005l
#define DEF_XB   0.000005l
#define DEF_YA  -0.000005l
#define DEF_YB   0.000005l
#define DEF_Z    0.00001l
#define DEF_W    0.00002l
#define DEF_L    0.00002l
#define DEF_I0   7'000 // like a 5 mW laser pointer
#define DEF_ATY  0.0625l/(1024.0l*1024.0l)
#define DEF_RTY  0.0625l/(1024.0l*1024.0l)
#define DEF_PTY  1/(1024.0l*1024.0l*1024.0l)
#define DEF_MDY  (uint64_t) 52
#define DEF_ATX  0.0625l/(1024.0l*1024.0l)
#define DEF_RTX  0.0625l/(1024.0l*1024.0l)
#define DEF_PTX  1/(1024.0l*1024.0l*1024.0l)
#define DEF_MDX  (uint64_t) 52

typedef union float_val {
    float f;
    double d;
    long double x;
} flt;

template <typename T>
struct sim_vals {
    uint64_t nx; // number of pixels along x-axis of detector
    uint64_t ny; // number of pixels along y-axis of detector
    T lam; // wavelength of light (m)
    T xa; // lower x-limit of aperture
    T xb; // upper x-limit of aperture
    T ya; // lower y-limit of aperture
    T yb; // upper y-limit of aperture
    T z; // distance to detector
    T w; // width of detector (along x-axis)
    T l; // length of detector (along y-axis)
    T I0; // intensity of incident light on aperture (W/m^2)
    T abstol_y;
    T reltol_y;
    T ptol_y;
    T mdepth_y;
    T abstol_x;
    T reltol_x;
    T ptol_x;
    T mdepth_x;
    uint64_t threads;
    uint64_t ptime;
    const char *bmp_path;
};

template <typename T, bool verbose> requires (std::is_floating_point_v<T>)
int start_sim(gtd::parser &parser) {
    sim_vals<T> vals;
    vals.nx = parser.get_arg("--nx", (uint64_t) DEF_NX);
    if (!vals.nx) {
        fprintf(stderr, "Error: number of pixels in detector along x-axis cannot be zero.\n");
        return 1;
    }
    vals.ny = parser.get_arg("--ny", (uint64_t) DEF_NY);
    if (!vals.ny) {
        fprintf(stderr, "Error: number of pixels in detector along y-axis cannot be zero.\n");
        return 1;
    }
    vals.lam = parser.get_arg("--lam", (T) DEF_LAM);
    if (vals.lam <= 0) {
        fprintf(stderr, "Error: wavelength of light must be positive.\n");
        return 1;
    }
    vals.xa = parser.get_arg("--xa", (T) DEF_XA);
    vals.xb = parser.get_arg("--xb", (T) DEF_XB);
    if (vals.xa >= vals.xb) {
        fprintf(stderr, "Error: lower x-limit of aperture must be less than upper x-limit of aperture.\n");
        return 1;
    }
    vals.ya = parser.get_arg("--ya", (T) DEF_YA);
    vals.yb = parser.get_arg("--yb", (T) DEF_YB);
    if (vals.ya >= vals.yb) {
        fprintf(stderr, "Error: lower y-limit of aperture must be less than upper y-limit of aperture.\n");
        return 1;
    }
    vals.z = parser.get_arg("--z", (T) DEF_Z);
    if (vals.z < 0) {
        fprintf(stderr, "Error: distance to detector along z-axis must be non-negative.\n");
        return 1;
    }
    vals.w = parser.get_arg("--w", (T) DEF_W);
    if (vals.w <= 0) {
        fprintf(stderr, "Error: width of detector (along x-axis) must be positive.\n");
        return 1;
    }
    vals.l = parser.get_arg("--l", (T) DEF_L);
    if (vals.l <= 0) {
        fprintf(stderr, "Error: length of detector (along y-axis) must be positive.\n");
        return 1;
    }
    vals.I0 = parser.get_arg("--I0", (T) DEF_I0);
    if (vals.I0 <= 0) {
        fprintf(stderr, "Error: intensity of incident light must be positive.\n");
        return 1;
    }
    vals.abstol_y = parser.get_arg("--abstol_y", (T) DEF_ATY);
    if (vals.abstol_y <= 0) {
        fprintf(stderr, "Error: absolute y-tolerance must be positive.\n");
        return 1;
    }
    vals.reltol_y = parser.get_arg("--reltol_y", (T) DEF_RTY);
    if (vals.reltol_y <= 0) {
        fprintf(stderr, "Error: relative y-tolerance must be positive.\n");
        return 1;
    }
    /* vals.ptol_y   = parser.get_arg("--ptol_y", (T) DEF_PTY);
    if (vals.ptol_y <= 0) {
        fprintf(stderr, "Error: periodic y-tolerance must be positive.\n");
        return 1;
    } */
    vals.abstol_x = parser.get_arg("--abstol_x", (T) DEF_ATX);
    if (vals.abstol_x <= 0) {
        fprintf(stderr, "Error: absolute x-tolerance must be positive.\n");
        return 1;
    }
    vals.reltol_x = parser.get_arg("--reltol_x", (T) DEF_RTX);
    if (vals.reltol_x <= 0) {
        fprintf(stderr, "Error: relative x-tolerance must be positive.\n");
        return 1;
    }
    /* vals.ptol_x   = parser.get_arg("--ptol_x", (T) 1/(1024.0l*1024.0l*1024.0l));
    if (vals.ptol_x <= 0) {
        fprintf(stderr, "Error: periodic x-tolerance must be positive.\n");
        return 1;
    } */
    diff::colourmap<T> cmap; // grayscale by default
    const char *cmap_str = parser.get_arg("--cmap");
    if (cmap_str) {
        if (auto it = diff::cmaps<T>::all_cmaps.find(cmap_str); it != diff::cmaps<T>::all_cmaps.end())
            cmap = it->second;
        else
            cmap.from_cmap(cmap_str);
    }
    // std::cout << cmap << std::endl;
    // bool prog = parser.get_arg("-p", true);
    vals.ptol_y   = parser.get_arg("--ptol_y", (T) DEF_PTY);
    vals.ptol_x   = parser.get_arg("--ptol_x", (T) DEF_PTX);
    vals.mdepth_y = parser.get_arg("--mdepth_y", DEF_MDY);
    vals.mdepth_x = parser.get_arg("--mdepth_x", DEF_MDX);
    vals.threads = parser.get_arg("--threads", (uint64_t) 0);
    vals.ptime = parser.get_arg("--ptime", (uint64_t) 1);
    vals.bmp_path = parser.get_arg("-o");
    std::cout << vals.ptol_y << ", " << vals.ptol_x << std::endl;
    if (!parser.empty()) {
        fprintf(stderr, "Error: unrecognised arguments have been passed:\n");
        for (const std::pair<int, std::string> &arg : parser)
            fprintf(stderr, "\t\"%s\"\n", arg.second.c_str());
        return 1;
    }
    if constexpr (verbose) {
        if constexpr (std::same_as<T, long double>)
            printf("Floating-point type: \"long double\"\n");
        else if (std::same_as<T, double>)
            printf("Floating-point type: \"double\"\n");
        else
            printf("Floating-point type: \"float\"\n");
    }
    vals.threads = vals.threads ? vals.threads : std::thread::hardware_concurrency();
    auto gfunc = [&vals](const T&){
        return vals.xa;
    };
    auto hfunc = [&vals](const T&){
        return vals.xb;
    };
    using G = decltype(gfunc);
    using H = decltype(hfunc);
    if constexpr (verbose)
        std::cout << "Starting simulation with the following parameters:"
                     "\n\tAperture lower x-limit = " << vals.xa       << " m"
                     "\n\tAperture upper x-limit = " << vals.xb       << " m"
                     "\n\tAperture lower y-limit = " << vals.ya       << " m"
                     "\n\tAperture upper y-limit = " << vals.yb       << " m"
                     "\n\tWavelength of light    = " << vals.lam      << " m"
                     "\n\tDistance to detector   = " << vals.z        << " m"
                     "\n\tWidth of detector      = " << vals.w        << " m"
                     "\n\tLength of detector     = " << vals.l        << " m"
                     "\n\tDetector x-resolution  = " << vals.nx       <<
                     "\n\tDetector y-resolution  = " << vals.ny       <<
                     "\n\tLight intensity        = " << vals.I0       << " W/m^2"
                     "\n\tAbsolute y-tolerance   = " << vals.abstol_y <<
                     "\n\tRelative y-tolerance   = " << vals.reltol_y <<
                     "\n\tPeriodic y-tolerance   = " << vals.ptol_y   <<
                     "\n\tAbsolute x-tolerance   = " << vals.abstol_x <<
                     "\n\tRelative x-tolerance   = " << vals.reltol_x <<
                     "\n\tPeriodic x-tolerance   = " << vals.ptol_x   <<
                     "\n\tMax. y-recursion-depth = " << vals.mdepth_y <<
                     "\n\tMax. x-recursion-depth = " << vals.mdepth_x <<
                     "\n\tNumber of threads      = " << vals.threads  <<
                     "\n\tProgress time delay    = " << vals.ptime    << " s\n";
    diff::aperture<T, G, H> ap{vals.ya, vals.yb, gfunc, hfunc};
    diff::diffimg<T, G, H> sim{vals.lam, ap, vals.z, vals.w, vals.l, vals.I0, vals.nx, vals.ny};
    sim.diffract(vals.abstol_y,
                 vals.reltol_y,
                 vals.ptol_y,
                 vals.mdepth_y,
                 vals.abstol_x,
                 vals.reltol_x,
                 vals.ptol_x,
                 vals.mdepth_x,
                 vals.threads,
                 vals.ptime);
    if constexpr (verbose)
        printf("Simulation completed, generating BMP...\n");
    off_t bmp_size = sim.gen_bmp(cmap, vals.bmp_path);
    if constexpr (verbose)
        printf("BMP generated with a size of %llu bytes.\n", (unsigned long long) bmp_size);
    return 0;
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *flt_type = parser.get_arg("-f");
    bool verbose = parser.get_arg("-v", false);
    if (verbose) {
        if (!flt_type || gtd::str_eq(flt_type, "Lf"))
            return start_sim<long double, true>(parser);
        if (gtd::str_eq(flt_type, "lf"))
            return start_sim<double, true>(parser);
#ifndef __CUDACC__
        if (gtd::str_eq(flt_type, "f"))
            return start_sim<float, true>(parser);
#endif
    } else {
        if (!flt_type || gtd::str_eq(flt_type, "Lf"))
            return start_sim<long double, false>(parser);
        if (gtd::str_eq(flt_type, "lf"))
            return start_sim<double, false>(parser);
#ifndef __CUDACC__
        if (gtd::str_eq(flt_type, "f"))
            return start_sim<float, false>(parser);
#endif
    }
    fprintf(stderr, "Error: invalid argument for \"-f\" flag, \"%s\".\n", flt_type);
    return 1;
}
