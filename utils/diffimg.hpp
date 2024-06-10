#ifndef DIFFIMG_HPP
#define DIFFIMG_HPP

#include "diffsim.hpp"
#include "../glib/nbod/gregbmp.hpp"
#include <map>

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
            typename std::map<T, colour<T>>::const_iterator lower = this->clrs.lower_bound(val);
            typename std::map<T, colour<T>>::const_iterator upper = this->clrs.upper_bound(val);
            T dtol = (val - lower->first)/(upper->first - lower->first);
            T dtou = 1 - dtol;
            return {dtol*lower->second.r + dtou*upper->second.r,
                    dtol*lower->second.g + dtou*upper->second.g,
                    dtol*lower->second.b + dtou*upper->second.b};
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
    template <gtd::numeric T, gtd::numeric C = long double>
    class diffimg : public diffsim<T> {
        gtd::bmp bmp{diffalloc<T>::nw, diffalloc<T>::nh};
        std::pair<T, T> minmax_val() const noexcept {
            T *ptr = diffalloc<T>::data;
            T minv = *ptr;
            T maxv = *ptr++;
            for (uint64_t i = 0; ++i < diffalloc<T>::np; ++ptr) {
                if (*ptr > maxv) {
                    maxv = *ptr;
                    continue;
                }
                if (*ptr < minv)
                    minv = *ptr;
            }
            return {minv, maxv};
        }
    public:
        using diffsim<T>::diffsim;
        void gen_bmp(const colourmap<C> &cmap = cmaps<T>::bgr, const char *path = nullptr, bool open = false) {
            uint64_t i = 0;
            uint64_t j = 0;
            T *ptr = diffalloc<T>::data;
            gtd::color color{};
            colour col{};
            auto [minv, maxv] = this->minmax_val();
            maxv -= minv;
            // T range = maxv - minv;
            while (i < diffalloc<T>::nw) {
                while (j < diffalloc<T>::nh) {
                    col = cmap((*ptr++ - minv)/maxv);
                    color.r = (unsigned char) std::round(col.r*255.0l);
                    color.g = (unsigned char) std::round(col.g*255.0l);
                    color.b = (unsigned char) std::round(col.b*255.0l);
                    bmp.set_pixel(i, j++, color);
                }
                ++i;
            }
            if (path)
                bmp.set_path(path);
            else
                bmp.set_path("name.bmp"); // make clever
            bmp.write();
            if (open)
                bmp.open_image();
        }
    };
}
#endif
