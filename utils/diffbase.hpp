#ifndef DIFFBASE_HPP
#define DIFFBASE_HPP

#include <iostream>
#include <cmath>
#include <cstdint>
#include <cinttypes>

namespace diff {
    template <typename T>
    concept numeric = requires (T a, T b) {
        {a + b} -> std::same_as<T>;
        {a + b} -> std::same_as<T>;
        {a * b} -> std::same_as<T>;
        {a / b} -> std::same_as<T>;
    };
    namespace misc {
        template <typename T = uint64_t> requires (std::is_integral<T>::value && std::is_unsigned<T>::value)
        T to_uint(const char *str) {
            if (!str || !*str)
                return -1;
            T res = 0;
            goto start;
            while (*str) {
                res *= 10;
                start:
                if (*str < 48 || *str > 57)
                    return -1;
                res += *str++ - 48;
            }
            return res;
        }
    }
    template <numeric T>
    T integrate_trap(T (*f)(const T&), const T &a, const T &b, uint64_t num = 1'000'000) {
        if (!f)
            return 0;
        T dx = (b - a)/num;
        uint64_t i = 0;
        uint64_t j = 1;
        T res = 0;
        while (i < num)
            res += (dx/2)*(f(a + i++*dx) + f(a + j++*dx));
        return res;
    }
}
#endif
