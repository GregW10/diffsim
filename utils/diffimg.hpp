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
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    struct colour {
        T r, g, b;
    };
#pragma pack(push, 1)
    struct bmp_col {
        unsigned char b, g, r;
    };
#pragma pack(pop)
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    class colours {
    public:
        constexpr static colour<T> white = {1.0l, 1.0l, 1.0l};
        constexpr static colour<T> black = {0.0l, 0.0l, 0.0l};
        constexpr static colour<T> blue  = {0.0l, 0.0l, 1.0l};
        constexpr static colour<T> green = {0.0l, 1.0l, 0.0l};
        constexpr static colour<T> red   = {1.0l, 0.0l, 0.0l};
    };
    template <typename T> requires (std::is_floating_point_v<T>)
    class colourmap {
        std::map<T, colour<T>> clrs{};
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
        colourmap(const colourmap<T> &other) : clrs{other.clrs} {}
        colourmap(colourmap<T> &&other) : clrs{std::move(other.clrs)} {}
        colour<T> operator()(const T &val) const {
            if (val < 0 || val > 1)
                throw std::invalid_argument{"Error: floating point value must be within [0,1].\n"};
            if (val == 0)
                return this->clrs.begin()->second;
            if (val == 1)
                return (--this->clrs.end())->second;
            // if (val > 0.5)
            //     std::cout << "val: " << val << std::endl;
            typename std::map<T, colour<T>>::const_iterator lower = this->clrs.lower_bound(val); // ge val
            typename std::map<T, colour<T>>::const_iterator upper = this->clrs.upper_bound(val); // gt val
            if (val < lower->first)
                --lower;
            // std::cout << "lower: " << lower->first << ", upper: " << upper->first << std::endl;
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
    };
    template <typename T = long double> requires (std::is_floating_point_v<T>)
    class cmaps {
    public:
        static inline const colourmap<T> grayscale{};
        static inline const colourmap<T> bgr = {{0.0l, colours<T>::blue},
                                                {0.5l, colours<T>::green},
                                                {1.0l, colours<T>::red}};
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
    template <gtd::numeric T>
    class diffimg : public diffsim<T> {
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
                gtd::write_all(fd, &(diffsim<T>::lambda), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::ap.ya), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::ap.yb), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::zdist), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::xdttr), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::ydttr), sizeof(T)) != sizeof(T) ||
                gtd::write_all(fd, &(diffsim<T>::E0), sizeof(T)) != sizeof(T))
                throw std::ios_base::failure{"Error: could not write .bmp metadata.\n"};
        }
    public:
        using diffsim<T>::diffsim;
        template <gtd::numeric C = T>
        off_t gen_bmp(const colourmap<C> &cmap = cmaps<T>::grayscale, const char *path = nullptr) {
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
                *oss << "diffpat_lam" << diffsim<T>::lambda << "m_ya"
                     << diffsim<T>::ap.ya << "m_yb" << diffsim<T>::ap.yb << "m_zd" << diffsim<T>::zdist
                     << "_xdl" << diffsim<T>::xdttr << "m_ydl" << diffsim<T>::ydttr << "m_E0"
                     << diffsim<T>::E0 << "Vpm.bmp";
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
