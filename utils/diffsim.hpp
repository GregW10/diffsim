#ifndef DIFFSIM_HPP
#define DIFFSIM_HPP

#include "diffalloc.hpp"
#include "../glib/misc/gregcomplex.hpp"
#include <cfloat>
#include <sstream>

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
    class invalid_dffr_format : public std::ios_base::failure {
    public:
        invalid_dffr_format() : std::ios_base::failure{"Error: invalid .dffr file format.\n"} {}
        explicit invalid_dffr_format(const char *msg) : std::ios_base::failure{msg} {}
    };
    class invalid_aperture_params : public std::logic_error {
    public:
        invalid_aperture_params() : std::logic_error{"Error: invalid aperture parameters.\n"} {}
        explicit invalid_aperture_params(const char *msg) : std::logic_error{msg} {}
    };
    template <gtd::numeric T>
    using xbfunc = T (*)(const T&);
#ifndef __CUDACC__
    template <gtd::numeric T = long double>//, gtd::callret<T> G = xbfunc<T>, gtd::callret<T> H = xbfunc<T>>
    requires (std::is_floating_point_v<T>)
#else
    template <gtd::numeric T = double>//, gtd::callret<T> G = xbfunc<T>, gtd::callret<T> H = xbfunc<T>>
            requires (std::is_floating_point_v<T> && !std::same_as<T, long double>) // no long double in GPU
#endif
    class diffsim;
#pragma pack(push, 1)
    template <typename T> requires (std::is_floating_point_v<T>)
    struct dffr_info {
        char hdr[4] = {'D', 'F', 'F', 'R'};
        uint32_t sizeT = (uint32_t) sizeof(T);
        uint32_t mant_dig = std::same_as<T, float> ? FLT_MANT_DIG : (std::same_as<T,double>?DBL_MANT_DIG:LDBL_MANT_DIG);
        T lam;
        T zd;
        T w;
        T l;
        T I_0;
        uint64_t nx;
        uint64_t ny;
    };
#pragma pack(pop)
    template <gtd::numeric T>//, /* gtd::callret<T> F, */ bool rec_calls = false>
    class functor {
        // F func;
    public:
        // explicit functor(const F &f) : func{f} {}
        // inline T operator()(const T &val) const {
        //     return this->func(val);
        // }
        virtual inline T operator()(const T&) const = 0;
    };
    /* template <gtd::numeric T, gtd::callret<T> F>
    class functor<T, F, true> {
        // F func;
        uint64_t ncalls{};
    public:
        explicit functor(const F &f) : func{f} {}
        inline T operator()(const T &val) {
            ++ncalls;
            return this->func(val);
        }
        uint64_t num_calls() const noexcept {
            return this->ncalls;
        }
        uint64_t reset() const noexcept {
            uint64_t num = this->ncalls;
            this->ncalls = 0;
            return num;
        }
    }; */
    template <gtd::numeric T>//, gtd::callret<T> G, gtd::callret<T> H>
    class aperture {
    protected:
        const uint32_t id;
        T ya;
        T yb;
        // G gfunc;
        // H hfunc;
    public:
        aperture(uint32_t ID, T _ya, T _yb) : id{ID}, ya{_ya}, yb{_yb} {} //, gfunc{_gfunc}, hfunc{_hfunc}{}
        virtual void write_ap_info(int fd) const = 0;
        virtual void read_ap_info(int fd) = 0;
        virtual void gen_fpath(const dffr_info<T> &inf, const char *suffix, std::string *out) const = 0;
        /* G gf() const noexcept {
            return this->gfunc;
        }
        H hf() const noexcept {
            return this->hfunc;
        } */
        virtual const functor<T> &gfunc() const noexcept = 0;
        virtual const functor<T> &hfunc() const noexcept = 0;
        virtual T gfunc(const T&) const = 0;
        virtual T hfunc(const T&) const = 0;
        friend class diffsim<T>;//, G, H>;
    };
    template <gtd::numeric T>//, gtd::callret<T> G, gtd::callret<T> H>
    class rectangle : public aperture<T> {//, T (*)(const T&), T (*)(const T&)> {
        // using GH = T (*)(const T&);
        T xa;
        T xb;
        class rc_functor : public functor<T> {
            T val;
        public:
            explicit rc_functor(const T &v) : val{v} {}
            inline T operator()(const T&) const {
                return this->val;
            }
        };
        rc_functor _gfunc{xa};
        rc_functor _hfunc{xb};
    public:
        static constexpr uint32_t rc_id = 0;
        // using aperture<T, G, H>::aperture;
        rectangle(T _xa, T _xb, T _ya, T _yb) : aperture<T>{rc_id, _ya, _yb}, xa{_xa}, xb{_xb} {}
        explicit rectangle(int fd, bool skip_id = false) {
            this->read_ap_info(fd, skip_id);
        }
        const functor<T> &gfunc() const noexcept {
            return this->_gfunc;
        }
        const functor<T> &hfunc() const noexcept {
            return this->_hfunc;
        }
        T gfunc(const T &val) const {
            return this->_gfunc(val);
        }
        T hfunc(const T &val) const {
            return this->_hfunc(val);
        }
        consteval static uint64_t ap_info_nb() noexcept {
            return sizeof(uint32_t) + 4*sizeof(T);
        }
        void write_ap_info(int fd) const {
            if (gtd::write_all(fd, &(aperture<T>::id), sizeof(uint32_t)) != sizeof(uint32_t))
                throw std::ios_base::failure{"Error: could not write aperture ID to .dffr file.\n"};
            // static constexpr uint32_t sizeT = (uint32_t) sizeof(T);
            T vals[4] = {this->xa, this->xb, this->ya, this->yb};
            if (gtd::write_all(fd, vals, 4*sizeof(T)) != 4*sizeof(T))
                throw std::ios_base::failure{"Error: could not write aperture limits to .dffr file.\n"};
        }
        void read_ap_info(int fd, bool skip_id = false) {
            if (!skip_id) {
                uint32_t _id;
                if (gtd::read_all(fd, &_id, sizeof(uint32_t)) != sizeof(uint32_t))
                    throw std::ios_base::failure{"Error: could not read aperture ID from .dffr file.\n"};
                if (_id != rc_id)
                    throw invalid_aperture_params{"Error: aperture ID in file does not match the ID of a rectangular "
                                                  "aperture.\n"};
            }
            // static constexpr uint32_t sizeT = (uint32_t) sizeof(T);
            T vals[4];
            if (gtd::read_all(fd, vals, 4*sizeof(T)) != 4*sizeof(T))
                throw std::ios_base::failure{"Error: could not read aperture limits from .dffr file.\n"};
            this->xa = vals[0];
            this->xb = vals[1];
            this->ya = vals[2];
            this->yb = vals[3];
        }
        void gen_fpath(const dffr_info<T> &inf, const char *suffix, std::string *out) const {
            if (!out)
                return;
            std::ostringstream oss{};
            oss << "diffpat_lam" << inf.lam << "m_xa" << this->xa << "m_xb" << this->xb << "m_ya"
                << this->ya << "m_yb" << this->yb << "m_zd" << inf.zd << "m_xdl" << inf.w << "m_ydl" << inf.l << "m_E0"
                << inf.I_0 << "Wpm2" << (suffix ? suffix : "");
            *out = oss.rdbuf()->view().data();
        }
        friend class diffsim<T>;
    };
    template <gtd::numeric T>
    struct vector {
        T x, y, z;
    };
#ifndef __CUDACC__
    template <gtd::numeric T>//, gtd::callret<T> G, gtd::callret<T> H>
            requires (std::is_floating_point_v<T>)
#else
    template <gtd::numeric T = double>//, gtd::callret<T> G = xbfunc<T>, gtd::callret<T> H = xbfunc<T>>
            requires (std::is_floating_point_v<T> && !std::same_as<T, long double>) // no long double in GPU
#endif
    class diffsim : public diffalloc<T> {
    protected:
        T lambda; // wavelength of light (m)
        const aperture<T> *ap; // aperture (slit)
        T zdist; // distance along z-axis to detector (m)
        T xdttr; // width of detector (along x) (m)
        T ydttr; // length of detector (along y) (m)
        T k; // wavenumber of light (rad m^-1)
        T I0; // intensity of light incident on aperture (W/m^2) - only stored to avoid recalculation
        T E0; // electric field amplitude of light incident on aperture (V/m)
        T pw; // pixel width (and height, since pixels are squares)
        T x0 = 0.5*(pw - xdttr);
        T y0 = 0.5*(pw - ydttr);
        T zdsq = zdist*zdist;
        gtd::complex<T> outside_factor = (gtd::complex<T>::m_imunit*zdist*E0)/lambda;
        static constexpr T eps0co2 = 0.5l*PERMITTIVITY*LIGHT_SPEED;
        bool ap_ff = false; // aperture from file - used to discern whether the aperture was loaded from a file or not
    private:
        std::atomic<uint64_t> counter{};
        void pix_to_pos(vector<T> *vec, uint64_t i, uint64_t j) {
            vec->x = x0 + i*pw;
            vec->y = y0 + j*pw;
        }
        void diffract_thread(T abstol_y,
                             T reltol_y,
                             T ptol_y,
                             uint64_t mdepth_y,
                             T abstol_x,
                             T reltol_x,
                             T ptol_x,
                             uint64_t mdepth_x) {
            coord c;
            vector<T> pos;
            pos.z = zdist; // is constant
            uint64_t offset = counter++;
            auto integrand = [&pos, this](const T &ap_x, const T &ap_y){
                T xdiff = ap_x - pos.x;
                T ydiff = ap_y - pos.y;
                T rsq = xdiff*xdiff + ydiff*ydiff + this->zdsq;
                T kr = this->k*std::sqrt(rsq);
                return ((std::cos(kr) + gtd::complex<T>::imunit*std::sin(kr))/
                       (rsq))*(1 + gtd::complex<T>::imunit/kr); // eliminate sine
            };
            while (offset < diffalloc<T>::np) {
                c = diffalloc<T>::pix_coords(offset); // I will find a more optimised way of doing this
                this->pix_to_pos(&pos, c.x, c.y);
                *(diffalloc<T>::data + offset) =
                        E0_to_intensity(diff::simpdblquad<T, gtd::complex<T>, decltype(integrand),
                                decltype(ap->gfunc()), decltype(ap->hfunc())>
                                         (integrand, ap->ya, ap->yb, ap->gfunc(), ap->hfunc(),
                                          abstol_y, reltol_y, ptol_y, &mdepth_y, abstol_x, reltol_x, ptol_x,
                                          &mdepth_x)*this->outside_factor); // get intensity
                offset = counter++;
            }
        }
        void prog_thread(unsigned ptime) {
            char _c = isatty(STDOUT_FILENO) ? '\r' : '\n';
            uint64_t offset = counter.load();
            while (offset < diffalloc<T>::np) {
                printf("%.2Lf%% complete%c", (((long double) offset)/diffalloc<T>::np)*100.0l, _c);
                fflush(stdout);
                sleep(ptime);
                offset = counter.load();
            }
            printf("100.00%% complete.\n");
        }
    public:
        static constexpr uint64_t def_mdepth = 52;
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
        ap{&aperture},
        zdist{dttr_dist},
        xdttr{dttr_width},
        ydttr{dttr_length},
        k{(T) (2*PI/wavelength)},
        pw{((T) dttr_width)/dttr_width_px},
        I0{incident_light_intensity},
        E0{intensity_to_E0(incident_light_intensity)}
        {
            if (wavelength <= 0 || aperture.ya >= aperture.yb/*||!aperture.gfunc||!aperture.hfunc*/|| dttr_dist < 0 ||
                dttr_width < 0 || dttr_length < 0)
                throw invalid_diffparams{};
        }
        void diffract(T abstol_y = 0.0625l/(1024.0l*1024.0l),
                      T reltol_y = 0.0625l/(1024.0l*1024.0l),
                      T ptol_y = 1/(1024.0l*1024.0l*1024.0l),
                      uint64_t mdepth_y = def_mdepth,
                      T abstol_x = 0.0625l/(1024.0l*1024.0l),
                      T reltol_x = 0.0625l/(1024.0l*1024.0l),
                      T ptol_x = 1/(1024.0l*1024.0l*1024.0l),
                      uint64_t mdepth_x = def_mdepth,
                      unsigned int num_threads = 0,
                      unsigned ptime = 1) {
            // be able to return max_depth overall
            static const unsigned int numt = std::thread::hardware_concurrency();
            std::vector<std::thread> threads;
            unsigned int _i = num_threads ? num_threads : numt;
            threads.reserve(_i);
            while (_i --> 0)
                threads.emplace_back(&diffsim<T>::diffract_thread, this, abstol_y, reltol_y, ptol_y, mdepth_y,
                                     abstol_x, reltol_x, ptol_x, mdepth_x);
            if (ptime) {
                std::thread{&diffsim<T>::prog_thread, this, ptime}.join();
                // progt.join();
            }
            for (std::thread &t : threads)
                t.join();
            counter.store(0);
        }
        off_t to_dffr(std::string &path) {
            static_assert(sizeof(T) <= ((uint32_t) -1), "\"sizeof(T)\" too large.\n"); // would never happen
            static_assert(LDBL_MANT_DIG <= ((uint32_t) -1), "\"LDBL_MANT_DIG\" too large.\n"); // would never happen
            dffr_info<T> info;
            info.lam = this->lambda;
            info.zd = this->zdist;
            info.w = this->xdttr;
            info.l = this->ydttr;
            info.I_0 = this->I0;
            info.nx = diffalloc<T>::nw;
            info.ny = diffalloc<T>::nh;
            if (path.empty())
                this->ap->gen_fpath(info, ".dffr", &path);
            int fd = open(path.c_str(), O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR);
            if (fd == -1)
                throw std::ios_base::failure{"Error: could not open .dffr file.\n"};
            if (gtd::write_all(fd, &info, sizeof(dffr_info<T>)) != sizeof(dffr_info<T>))
                throw std::ios_base::failure{"Error: could not write .dffr file information.\n"};
            this->ap->write_ap_info(fd);
            if (gtd::write_all(fd, this->data, this->nb) != this->nb)
                throw std::ios_base::failure{"Error: could not write 2D intensity data array to .dffr file.\n"};
            off_t _end = lseek(fd, 0, SEEK_CUR);
            if (close(fd) == -1)
                throw std::ios_base::failure{"Error: could not close .dffr file.\n"};
            return _end;
        }
        uint64_t from_dffr(const char *path) {
            static_assert(sizeof(T) <= ((uint32_t) -1), "\"sizeof(T)\" too large.\n"); // would never happen
            static_assert(LDBL_MANT_DIG <= ((uint32_t) -1), "\"LDBL_MANT_DIG\" too large.\n"); // would never happen
            struct stat buff{};
            if (stat(path, &buff) == -1)
                throw std::ios_base::failure{"Error: could not obtain file information.\n"};
            if (!S_ISREG(buff.st_mode))
                throw std::ios_base::failure{"Error: file is not a regular file.\n"};
            if (buff.st_size < sizeof(dffr_info<T>) + sizeof(uint32_t))
                throw invalid_dffr_format{"Error: .dffr file size too small.\n"};
            int fd = open(path, O_RDONLY);
            if (fd == -1)
                throw std::ios_base::failure{"Error: could not open .dffr file.\n"};
            dffr_info<T> info;
            if (gtd::read_all(fd, &info, sizeof(dffr_info<T>)) != sizeof(dffr_info<T>))
                throw std::ios_base::failure{"Error: could not read .dffr file information.\n"};
            if (info.hdr[0] != 'D' || info.hdr[1] != 'F' || info.hdr[2] != 'F' || info.hdr[3] != 'R')
                throw invalid_dffr_format{"Error: invalid .dffr file header.\n"};
            if (info.sizeT != sizeof(T))
                throw invalid_dffr_format{"Error: reported \"sizeof(T)\" does not match actual \"sizeof(T)\".\n"};
            if constexpr (std::same_as<T, long double>) {
                if (info.mant_dig != LDBL_MANT_DIG)
                    throw invalid_dffr_format{"Error: reported \"LDBL_MANT_DIG\" does not match actual "
                                              "\"LDBL_MANT_DIG\".\n"};
            } else if constexpr (std::same_as<T, double>) {
                if (info.mant_dig != DBL_MANT_DIG)
                    throw invalid_dffr_format{"Error: reported \"DBL_MANT_DIG\" does not match actual "
                                              "\"DBL_MANT_DIG\".\n"};
            } else {
                if (info.mant_dig != FLT_MANT_DIG)
                    throw invalid_dffr_format{"Error: reported \"FLT_MANT_DIG\" does not match actual "
                                              "\"FLT_MANT_DIG\".\n"};
            }
            if (info.lam <= 0)
                throw invalid_dffr_format{"Error: wavelength is not positive.\n"};
            if (info.zd < 0)
                throw invalid_dffr_format{"Error: z-distance is negative.\n"};
            if (info.w < 0)
                throw invalid_dffr_format{"Error: detector width is negative.\n"};
            if (info.l < 0)
                throw invalid_dffr_format{"Error: detector length is negative.\n"};
            if (info.I_0 < 0)
                throw invalid_dffr_format{"Error: incident light intensity is negative.\n"};
            uint32_t ap_id;
            if (gtd::read_all(fd, &ap_id, sizeof(uint32_t)) != sizeof(uint32_t))
                throw std::ios_base::failure{"Error: could not read aperture ID from .dffr file.\n"};
            uint64_t psize;
            if (ap_id == rectangle<T>::rc_id) {
                if (buff.st_size < (psize = sizeof(dffr_info<T>) + rectangle<T>::ap_info_nb()))
                    throw invalid_dffr_format{"Error: insufficient file size for rectangular aperture.\n"};
                this->ap = new rectangle<T>{fd, true}; // skip ID, since already read
                ap_ff = true;
            } else {
                throw invalid_aperture_params{"Error: unknown aperture ID encountered.\n"};
            }
            // this->ap->read_ap_info(fd, true);
            diffalloc<T>::np = info.nx*info.ny;
            diffalloc<T>::nb = diffalloc<T>::np*sizeof(T);
            // psize += info.nx*info.ny*sizeof(T);
            if (buff.st_size != psize + diffalloc<T>::nb)
                throw invalid_dffr_format{"Error: invalid .dffr file size.\n"};
            diffalloc<T>::mapper.reset(diffalloc<T>::nb);
            diffalloc<T>::data = (T*) diffalloc<T>::mapper.get();
            if (gtd::read_all(fd, diffalloc<T>::data, diffalloc<T>::nb) != this->nb)
                throw std::ios_base::failure{"Error: could not read 2D intensity data array from .dffr file.\n"};
            // off_t _end = lseek(fd, 0, SEEK_CUR);
            if (close(fd) == -1)
                throw std::ios_base::failure{"Error: could not close .dffr file.\n"};
            diffalloc<T>::nw = info.nx;
            diffalloc<T>::nh = info.ny;
            this->lambda = info.lam;
            this->zdist = info.zd;
            this->xdttr = info.w;
            this->ydttr = info.l;
            this->k = (2*PI)/info.lam;
            this->I0 = info.I_0;
            this->E0 = intensity_to_E0(info.I_0);
            this->pw = info.w/info.nx;
            this->x0 = 0.5*(this->pw - info.w);
            this->y0 = 0.5*(this->pw - info.l);
            this->zdsq = info.zd*info.zd;
            this->outside_factor = (gtd::complex<T>::m_imunit*info.zd*this->E0)/info.lam;
            return psize;
        }
        ~diffsim() {
            if (ap_ff)
                delete this->ap;
            ap_ff = false; // overkill
        }
    };
}
#endif
