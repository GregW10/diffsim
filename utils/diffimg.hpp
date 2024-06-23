#ifndef DIFFIMG_HPP
#define DIFFIMG_HPP

#include "diffsim.hpp"
// #include "../glib/nbod/gregbmp.hpp"
#include <map>
#include <sstream>

namespace diff {
    class invalid_colour_format : public std::invalid_argument {
    public:
        invalid_colour_format() : std::invalid_argument{"Error: invalid colour format.\n"} {}
        explicit invalid_colour_format(const char *msg) : std::invalid_argument{msg} {}
    };
    class invalid_cmap_format : public std::ios_base::failure {
    public:
        invalid_cmap_format() : std::ios_base::failure{"Error: invalid .cmap file format.\n.\n"} {}
        explicit invalid_cmap_format(const char *msg) : std::ios_base::failure{msg} {}
    };
#pragma pack(push, 1)
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    struct colour {
        T r, g, b;
    };
    struct bmp_col {
        unsigned char b, g, r;
    };
#pragma pack(pop)
    template <typename T> requires (std::is_floating_point_v<T>)
    std::ostream &operator<<(std::ostream &os, const colour<T> &col) {
        return os << "(r=" << col.r << ",g=" << col.g << ",b=" << col.b << ')';
    }
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    class colours {
    public:
        constexpr static colour<T> white   = {1.0l, 1.0l, 1.0l};
        constexpr static colour<T> black   = {0.0l, 0.0l, 0.0l};
        constexpr static colour<T> blue    = {0.0l, 0.0l, 1.0l};
        constexpr static colour<T> green   = {0.0l, 1.0l, 0.0l};
        constexpr static colour<T> red     = {1.0l, 0.0l, 0.0l};
        constexpr static colour<T> magenta = {1.0l, 0.0l, 1.0l};
        constexpr static colour<T> cyan    = {0.0l, 1.0l, 1.0l};
        static inline const std::map<const char*, colour<T>> all_cols = {
                {"white", white},
                {"black", black},
                {"blue", blue},
                {"green", green},
                {"red", red},
                {"magenta", magenta},
                {"cyan", cyan}
        };
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class colourmap {
        std::map<T, colour<T>> clrs{};
#pragma pack(push, 1)
        struct cmap_el {
            T val;
            colour<T> col;
        };
#pragma pack(pop)
    public:
        colourmap() : clrs{{0, colours<T>::black}, {1, colours<T>::white}} {}
        colourmap(const std::initializer_list<std::pair<const T, colour<T>>> &list) /*: clrs{list}*/ {
            const typename std::initializer_list<std::pair<const T, colour<T>>>::size_type lsize = list.size();
            if (lsize < 2)
                throw std::invalid_argument{"Error: at least two colours must be passed.\n"};
            const std::pair<const T, colour<T>> *data = list.begin();
            if (data->first != 0 || (data + lsize - 1)->first != 1)
                throw invalid_colour_format{"Error: first and last colours in list must correspond to floating point "
                                            "values of 0 and 1, respectively.\n"};
            clrs = list;
        }
        colourmap(const std::pair<const T, colour<T>> *_clrs, uint64_t _size) /*: clrs{list}*/ {
            if (!_clrs)
                throw std::invalid_argument{"Error: \"nullptr\" cannot be passed as \"_clrs\" argument.\n"};
            if (_size < 2)
                throw std::invalid_argument{"Error: at least two colours must be passed.\n"};
            std::pair<const T, colour<T>> *_end = _clrs + _size;
            if (_clrs->first != 0 || (_end - 1)->first != 1)
                throw invalid_colour_format{"Error: first and last colours in list must correspond to floating point "
                                            "values of 0 and 1, respectively.\n"};
            clrs.insert(_clrs, _end);
        }
        explicit colourmap(const char *cmap_path) {
            this->from_cmap(cmap_path);
        }
        colourmap(const colourmap<T> &other) : clrs{other.clrs} {}
        colourmap(colourmap<T> &&other) : clrs{std::move(other.clrs)} {}
        off_t to_cmap(const char *path = "colourmap.cmap") const {
            if (!path || !*path)
                throw std::invalid_argument{"Error: a path MUST be passed.\n"};
            int fd = open(path, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR);
            if (fd == -1)
                throw std::ios_base::failure{"Error: could not open file.\n"};
            char buff[4] = {'C', 'M', 'A', 'P'};
            if (gtd::write_all(fd, buff, 4*sizeof(char)) != 4*sizeof(char))
                throw std::ios_base::failure{"Error: could not write header to .cmap file.\n"};
            uint64_t num_colours = clrs.size();
            if (gtd::write_all(fd, &num_colours, sizeof(uint64_t)) != sizeof(uint64_t))
                throw std::ios_base::failure{"Error: could not write the number of colours to .cmap file.\n"};
            static constexpr uint64_t sizeT = sizeof(T);
            if (gtd::write_all(fd, &sizeT, sizeof(uint64_t)) != sizeof(uint64_t))
                throw std::ios_base::failure{"Error: could not write sizeof(T) to .cmap file.\n"};
            cmap_el _el;
            for (const std::pair<T, colour<T>> &_p : this->clrs) {
                _el.val = _p.first;
                _el.col = _p.second;
                if (gtd::write_all(fd, &_el, sizeof(cmap_el)) != sizeof(cmap_el))
                    throw std::ios_base::failure{"Error: could not write colours array to .cmap file.\n"};
            }
            off_t pos = lseek(fd, 0, SEEK_CUR);
            if (close(fd) == -1)
                throw std::ios_base::failure{"Error: could not close file.\n"};
            return pos;
        }
        off_t from_cmap(const char *path) {
            if (!path || !*path)
                throw std::invalid_argument{"Error: a path MUST be passed.\n"};
            struct stat buffer{};
            if (stat(path, &buffer) == -1)
                throw std::ios_base::failure{"Error: could not obtain .cmap file information.\n"};
            if (!S_ISREG(buffer.st_mode))
                throw std::ios_base::failure{"Error: .cmap file is not a regular file.\n"};
            constexpr static uint64_t info_hdr_size = 4*sizeof(char) + 2*sizeof(uint64_t);
            if (buffer.st_size < info_hdr_size)
                throw invalid_cmap_format{"Error: insufficient .cmap file size.\n"};
            int fd = open(path, O_RDONLY);
            if (fd == -1)
                throw std::ios_base::failure{"Error: could not open file.\n"};
            char buff[5];
            buff[4] = 0;
            if (gtd::read_all(fd, buff, 4*sizeof(char)) != 4*sizeof(char))
                throw std::ios_base::failure{"Error: could not read header of .cmap file.\n"};
            if (!gtd::str_eq(buff, "CMAP"))
                throw invalid_cmap_format{"Error: invalid .cmap header.\n"};
            uint64_t num_colours;
            if (gtd::read_all(fd, &num_colours, sizeof(uint64_t)) != sizeof(uint64_t))
                throw std::ios_base::failure{"Error: could not read the number of colours from .cmap file.\n"};
            if (num_colours < 2)
                throw invalid_cmap_format{"Error: a colourmap cannot have less than two colours.\n"};
            uint64_t sizeT;
            if (gtd::read_all(fd, &sizeT, sizeof(uint64_t)) != sizeof(uint64_t))
                throw std::ios_base::failure{"Error: could not read sizeof(T) from .cmap file.\n"};
            if (sizeT != sizeof(T))
                throw invalid_cmap_format{"Error: sizeof(T) in .cmap file does not match actual sizeof(T).\n"};
            if (buffer.st_size != info_hdr_size + num_colours*sizeof(cmap_el))
                throw invalid_cmap_format{"Error: invalid .cmap file size.\n"};
            cmap_el _el;
            std::pair<T, colour<T>> _p;
            typename std::map<T, colour<T>>::iterator _end = this->clrs.end();
            while (num_colours --> 0) {
                if (gtd::read_all(fd, &_el, sizeof(cmap_el)) != sizeof(cmap_el))
                    throw std::ios_base::failure{"Error: could not write colours array to .cmap file.\n"};
                if (_el.val < 0 || _el.val > 1)
                    throw invalid_cmap_format{"Error: value in .cmap found outside [0,1] range.\n"};
                if (_el.col.r < 0 || _el.col.r > 1)
                    throw invalid_cmap_format{"Error: red value in .cmap found outside [0,1] range.\n"};
                if (_el.col.g < 0 || _el.col.g > 1)
                    throw invalid_cmap_format{"Error: green value in .cmap found outside [0,1] range.\n"};
                if (_el.col.b < 0 || _el.col.b > 1)
                    throw invalid_cmap_format{"Error: blue value in .cmap found outside [0,1] range.\n"};
                _p.first = _el.val;
                _p.second = _el.col;
                this->clrs.insert(_end, _p);
            }
            off_t pos = lseek(fd, 0, SEEK_CUR);
            if (close(fd) == -1)
                throw std::ios_base::failure{"Error: could not close file.\n"};
            return pos;
        }
        colour<T> operator()(const T &val) const {
            if (val < 0 || val > 1)
                throw std::invalid_argument{"Error: floating point value must be within [0,1].\n"};
            if (val == 0)
                return this->clrs.begin()->second;
            if (val == 1)
                return (--this->clrs.end())->second;
            typename std::map<T, colour<T>>::const_iterator lower = this->clrs.lower_bound(val); // ge val
            typename std::map<T, colour<T>>::const_iterator upper = this->clrs.upper_bound(val); // gt val
            if (val < lower->first)
                --lower;
            T dtol = (val - lower->first)/(upper->first - lower->first);
            T dtou = 1 - dtol;
            return {dtou*lower->second.r + dtol*upper->second.r,
                    dtou*lower->second.g + dtol*upper->second.g,
                    dtou*lower->second.b + dtol*upper->second.b};
        }
        colourmap<T> &operator=(const colourmap<T> &other) {
            if (&other == this)
                return *this;
            this->clrs = other.clrs;
            return *this;
        }
        colourmap<T> &operator=(colourmap<T> &&other) {
            if (&other == this)
                return *this;
            this->clrs = std::move(other.clrs);
            return *this;
        }
        // template <typename U>
        friend std::ostream &operator<<(std::ostream &os, const colourmap<T> &_cmap) {
            uint64_t _num = _cmap.clrs.size();
            if (!_num) // impossible, but anyway...
                return os << "[]";
            typename std::map<T, colour<T>>::const_iterator it = _cmap.clrs.cbegin();
            os << "[\n";
            goto start;
            while (--_num > 0) {
                os << ",\n";
                start:
                os << '\t' << it->first << ':' << it++->second;
            }
            return os << "\n]";
        }
    };
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    class cmaps {
    public:
        static inline const colourmap<T> grayscale{};
        static inline const colourmap<T> bgr = {{0.0l, colours<T>::blue},
                                                {0.5l, colours<T>::green},
                                                {1.0l, colours<T>::red}};
        static inline const std::map<std::string, colourmap<T>> all_cmaps = {
                {"grayscale", grayscale},
                {"bgr", bgr}
        };
    };
#pragma pack(push, 1)
    struct bmp_header {
        char hdr[2] = {'B', 'M'};
        uint32_t fsize{};
        uint16_t res1[2]{};
        uint32_t arr_offset = 54;
        uint32_t hsize = 40;
        uint32_t width{};
        uint32_t height{};
        uint16_t num_cpanes = 1;
        uint16_t bpp = 24;
        uint32_t res2[6]{};
    };
#pragma pack(pop)
#ifndef __CUDACC__
    template <gtd::numeric T = long double, gtd::callret<T> G = xbfunc<T>, gtd::callret<T> H = xbfunc<T>>
#else
    template <gtd::numeric T = double, gtd::callret<T> G = xbfunc<T>, gtd::callret<T> H = xbfunc<T>>
#endif
    class diffimg : public diffsim<T, G, H> {
        // gtd::bmp bmp{diffalloc<T>::nw, diffalloc<T>::nh};
        // gtd::mmapper bmp{diffalloc<T>::nb};
        std::pair<T, T> minmax_vals() const noexcept {
            T *ptr = diffalloc<T>::data;
            // std::cout << "*ptr: " << *ptr << std::endl;
            T minv = *ptr;
            T maxv = *ptr++;
            for (uint64_t i = 0; ++i < diffalloc<T>::np; ++ptr) {
                // if (*ptr > 0.001)
                //     std::cout << "*ptr: " << *ptr << std::endl;
                if (*ptr > maxv) {
                    maxv = *ptr;
                    continue;
                }
                if (*ptr < minv)
                    minv = *ptr;
            }
            // std::cout << "MINV: " << minv << ", MAXV: " << maxv << std::endl;
            return {minv, maxv};
        }
        void add_bmp_metadata(int fd) {
            constexpr static uint64_t sizeT = sizeof(T);
            if (gtd::write_all(fd, &sizeT, sizeof(uint64_t)) != sizeof(uint64_t) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::lambda), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::ap.ya), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::ap.yb), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::zdist), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::xdttr), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::ydttr), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T, G, H>::E0), sizeof(T)) != sizeof(T))
                throw std::ios_base::failure{"Error: could not write .bmp metadata.\n"};
        }
    public:
        using diffsim<T, G, H>::diffsim;
        template <gtd::numeric C = T>
        off_t gen_bmp(const colourmap<C> &cmap = cmaps<T>::bgr, const char *path = nullptr) {
            uint64_t i = 0;
            uint64_t j;
            T *dptr = diffalloc<T>::data;
            // gtd::color color{};
            colour<T> col{};
            auto [minv, maxv] = this->minmax_vals();
            maxv -= minv;
            // T range = maxv - minv;
            uint64_t wb = diffalloc<T>::nw*3;
            uint8_t rem = wb % 4; // row must end on 4-byte boundary:
            uint8_t padding = rem ? 4 - rem : 0;
            uint64_t wpb = wb + padding;
            uint8_t pad[3]{};
            uint64_t bmp_arrsize = diffalloc<T>::np*3 + diffalloc<T>::nh*padding;
            gtd::mmapper bmp{bmp_arrsize};
            unsigned char *uptr = (unsigned char*) bmp.get();
            bmp_col *cptr;
            while (i < diffalloc<T>::nw) {
                cptr = (bmp_col*) uptr;
                j = 0;
                while (j < diffalloc<T>::nh) {
                    // std::cout << "*dptr: " << *dptr << ", minv: " << minv << ", maxv: " << maxv << std::endl;
                    // if ((*dptr - minv)/maxv > 0.5)
                    //     std::cout << "argument: " << (*dptr - minv)/maxv << std::endl;
                    // if (*dptr > 0.001)
                    //     std::cout << "here: " << *dptr << std::endl;
                    col = cmap((*dptr++ - minv)/maxv);
                    // std::cout << col.r << ", " << col.g << ", " << col.b << std::endl;
                    cptr->r = (unsigned char) std::round(col.r*255.0l);
                    cptr->g = (unsigned char) std::round(col.g*255.0l);
                    cptr->b = (unsigned char) std::round(col.b*255.0l);
                    // color.r = (unsigned char) std::round(col.r*255.0l);
                    // color.g = (unsigned char) std::round(col.g*255.0l);
                    // color.b = (unsigned char) std::round(col.b*255.0l);
                    // bmp.set_pixel(i, j++, color);
                    ++cptr;
                    ++j;
                }
                uptr += wpb;
                ++i;
            }
            // char *fpath = const_cast<char*>(path);
            std::ostringstream *oss{};
            if (!path) {
                /* // rem = 255;
                errno = 0;
                long lim = pathconf(".", _PC_NAME_MAX);
                if (lim == -1 && errno)
                    throw std::ios_base::failure{"Error: could not obtain \"pathconf\" info.\n"};
                if (!lim)
                    lim = 1023;
                fpath = new char[lim + 1];
                snprintf(fpath, lim + 1, "difpat_lam"); */
                oss = new std::ostringstream{};
                *oss << "diffpat_lam" << diffsim<T, G, H>::lambda << "m_ya"
                     << diffsim<T, G, H>::ap.ya << "m_yb" << diffsim<T, G, H>::ap.yb << "m_zd"
                     << diffsim<T, G, H>::zdist
                     << "_xdl" << diffsim<T, G, H>::xdttr << "m_ydl" << diffsim<T, G, H>::ydttr << "m_E0"
                     << diffsim<T, G, H>::E0 << "Vpm.bmp";
                path = oss->rdbuf()->view().data();
            }
            /* if (path)
                bmp.set_path(path);
            else
                bmp.set_path("name.bmp"); // make clever
            bmp.write();
            if (open)
                bmp.open_image(); */
            int fd = open(path, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
            if (oss)
                delete oss;
            if (fd == -1)
                throw std::ios_base::failure{"Error: could not open/create .bmp file.\n"};
            bmp_header header;
            header.fsize = sizeof(header) + bmp_arrsize;
            header.width = diffalloc<T>::nw;
            header.height = diffalloc<T>::nh;
            if (gtd::write_all(fd, &header, sizeof(header)) != sizeof(header))
                throw std::ios_base::failure{"Error: could not write .bmp header.\n"};
            if (gtd::write_all(fd, bmp.get(), bmp_arrsize) != bmp_arrsize)
                throw std::ios_base::failure{"Error: could not write .bmp colour array.\n"};
            this->add_bmp_metadata(fd);
            off_t pos = lseek(fd, 0, SEEK_CUR);
            if (close(fd) == -1)
                throw std::ios_base::failure{"Error: could not close .bmp file.\n"};
            return pos;
        }
    };
}
#endif
