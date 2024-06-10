#ifndef DIFFSIM_HPP
#define DIFFSIM_HPP

#include "diffalloc.hpp"
#include "../glib/misc/gregcomplex.hpp"

#ifndef __CUDACC__
#include <thread>
#include <atomic>
#endif

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
    protected:
        T lambda; // wavelength of light (m)
        aperture<T> ap; // aperture (slit)
        T zdist; // distance along z-axis to detector (m)
        T xdttr; // width of detector (along x) (m)
        T ydttr; // length of detector (along y) (m)
        T k; // wavenumber of light (rad m^-1)
        T E0; // electric field amplitude of light incident on aperture (V/m)
        T pw; // pixel width (and height, since pixels are squares)
        T x0 = 0.5*(pw - xdttr);
        T y0 = 0.5*(pw - ydttr);
        static constexpr T eps0co2 = 0.5l*PERMITTIVITY*LIGHT_SPEED;
    private:
        std::atomic<uint64_t> counter{};
        void pix_to_pos(vector<T> *vec, uint64_t i, uint64_t j) {
            vec->x = x0 + i*pw;
            vec->y = y0 + j*pw;
        }
        void diffract_thread(T abstol_y,
                             T reltol_y,
                             T ptol_y,
                             uint64_t *mdepth_y,
                             T abstol_x,
                             T reltol_x,
                             T ptol_x,
                             uint64_t *mdepth_x) {
            coord c;
            vector<T> pos;
            pos.z = zdist; // is constant
            T zdsq = zdist*zdist;
            uint64_t offset = counter++;
            static const gtd::complex<T> _i{0, 1};
            // gtd::complex<T> result;
            auto integrand = [&pos, this, &zdsq](const T &ap_x, const T &ap_y){
                T xdiff = ap_x - pos.x;
                T ydiff = ap_y - pos.y;
                T rsq = xdiff*xdiff + ydiff*ydiff + zdsq;
                T kr = this->k*std::sqrt(rsq);
                return ((std::cos(kr) + _i*std::sin(kr))/(rsq))*(1 + _i/kr); // eliminate sine
            };
            while (offset < diffalloc<T>::np) {
                c = diffalloc<T>::pix_coords(offset); // I will find a more optimised way of doing this
                this->pix_to_pos(&pos, c.x, c.y);
                *(diffalloc<T>::data + offset) = E0_to_intensity(diff::simpdblquad<T, gtd::complex<T>,
                    decltype(integrand)>(integrand, ap.ya, ap.yb, ap.gfunc,
                                         ap.hfunc, abstol_y,
                                         reltol_y, ptol_y, mdepth_y, abstol_x, reltol_x, ptol_x,
                                         mdepth_x)*(((-_i)*zdist*E0)/lambda)); // get intensity
                offset = counter++;
            }
        }
    public:
        static constexpr uint64_t def_mdepth = 16;
        diffsim() = delete;
        diffsim(const T &wavelength,
                const aperture<T> &aperture,
                const T &dttr_dist,
                const T &dttr_width,
                const T &dttr_length,
                const T &incident_light_intensity,
                uint64_t dttr_width_px = 2'000,
                uint64_t dttr_length_px = 2'000) :
        diffalloc<T>{dttr_width_px, dttr_length_px},
        lambda{wavelength},
        ap{aperture},
        zdist{dttr_dist},
        xdttr{dttr_width},
        ydttr{dttr_length},
        k{2*PI/wavelength},
        pw{((long double) dttr_width)/dttr_width_px},
        E0{intensity_to_E0(incident_light_intensity)}
        {
            if (wavelength <= 0 || aperture.ya >= aperture.yb || !aperture.gfunc || !aperture.hfunc || dttr_dist < 0 ||
                dttr_width < 0 || dttr_length < 0)
                throw invalid_diffparams{};
        }
        void diffract(T abstol_y = 0.0625l/(1024.0l*1024.0l*1024.0l),
                      T reltol_y = 0.0625l/(1024.0l*1024.0l*1024.0l),
                      T ptol_y = 1/(1024.0l*1024.0l*1024.0l),
                      uint64_t mdepth_y = def_mdepth,
                      T abstol_x = 0.0625l/(1024.0l*1024.0l*1024.0l),
                      T reltol_x = 0.0625l/(1024.0l*1024.0l*1024.0l),
                      T ptol_x = 1/(1024.0l*1024.0l*1024.0l),
                      uint64_t mdepth_x = def_mdepth,
                      unsigned int num_threads = 0) {
            static const unsigned int numt = std::thread::hardware_concurrency();
            std::vector<std::thread> threads;
            unsigned int i = num_threads ? num_threads : numt;
            while (i --> 0)
                threads.emplace_back(&diffsim<T>::diffract_thread, this, abstol_y, reltol_y, ptol_y, &mdepth_y,
                                     abstol_x, reltol_x, ptol_x, &mdepth_x);
            for (std::thread &t : threads)
                t.join();
            counter.store(0);
        }
    };
}
#endif
