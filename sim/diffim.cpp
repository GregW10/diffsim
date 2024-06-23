#include "../utils/diffimg.hpp"
#include "../glib/misc/gregparse.hpp"

#define DEF_NX 500
#define DEF_NY 500
#define DEF_LAM 0.000000450l // blue-ish purple light
#define DEF_XA -0.000001l
#define DEF_XB 0.000001l
#define DEF_YA -0.000001l
#define DEF_YB 0.000001l
#define DEF_Z 0.000001l
#define DEF_W 0.000002l
#define DEF_L 0.000002l

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
};


template <typename T> requires (std::is_floating_point_v<T>)
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
    diff::colourmap<T> cmap; // grayscale by default
    const char *cmap_str = parser.get_arg("--cmap");
    if (cmap_str) {
        if (auto it = diff::cmaps<T>::all_cmaps.find(cmap_str); it != diff::cmaps<T>::all_cmaps.end())
            cmap = it->second;
        else
            cmap.from_cmap(cmap_str);
    }
    std::cout << cmap << std::endl;
    if (!parser.empty()) {
        fprintf(stderr, "Error: unrecognised arguments have been passed:\n");
        for (const std::pair<int, std::string> &arg : parser)
            fprintf(stderr, "\t\"%s\"\n", arg.second.c_str());
        return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *flt_type = parser.get_arg("-f");
    if (!flt_type || gtd::str_eq(flt_type, "f"))
        return start_sim<float>(parser);
    if (gtd::str_eq(flt_type, "lf"))
        return start_sim<double>(parser);
#ifndef __CUDACC__
    if (gtd::str_eq(flt_type, "Lf"))
        return start_sim<long double>(parser);
#endif
    fprintf(stderr, "Error: invalid argument for \"-f\" flag, \"%s\".\n", flt_type);
    return 1;
}
