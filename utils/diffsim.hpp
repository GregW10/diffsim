#ifndef DIFFSIM_HPP
#define DIFFSIM_HPP

#include "diffalloc.hpp"

namespace diff {
    class invalid_diffparams : public std::logic_error {
    public:
        invalid_diffparams() : std::logic_error{"Error: invalid diffraction simulation parameters.\n"} {}
        explicit invalid_diffparams(const char *msg) : std::logic_error{msg} {}
    };
    template <gtd::numeric T>
    using xbfunc = T (*)(const T&);
    template <gtd::numeric T>
    struct aperture {
        T ya;
        T yb;
        xbfunc<T> gfunc;
        xbfunc<T> hfunc;
    };
    template <gtd::numeric T>
    struct vector {
        T x, y, z;
    };
    template <gtd::numeric T>
    class diffsim : public diffalloc<T> {
        T lambda; // wavelength of light
        aperture<T> ap; // aperture (slit)
        T zdist; // distance along z-axis to detector
        T xdttr; // width of detector (along x)
        T ydttr; // length of detector (along y)
        T k; // wavenumber of light
        vector<T> pix_to_pos(uint64_t i, uint64_t j) {

        }
    public:
        diffsim() = delete;
        diffsim(const T &wavelength,
                const aperture<T> &aperture,
                const T &dttr_dist,
                const T &dttr_width,
                const T &dttr_length,
                uint64_t dttr_width_px = 2'000,
                uint64_t dttr_length_px = 2'000) :
        diffalloc<T>{dttr_width_px, dttr_length_px},
        lambda{wavelength},
        ap{aperture},
        zdist{dttr_dist},
        xdttr{dttr_width},
        ydttr{dttr_length},
        k{2*PI/wavelength} {
            if (wavelength <= 0 || aperture.ya >= aperture.yb || !aperture.gfunc || !aperture.hfunc || dttr_dist < 0 ||
                dttr_width < 0 || dttr_length < 0)
                throw invalid_diffparams{};
        }

    };
}
#endif
