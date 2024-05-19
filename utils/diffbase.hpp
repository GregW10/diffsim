#ifndef DIFFBASE_HPP
#define DIFFBASE_HPP

#if !defined(__linux__) && !defined(__APPLE__)
#error OS not supported
#endif

#include <climits>

#if (CHAR_BIT != 8)
#error A byte has to have 8 bits.
#endif

#ifdef __CUDACC__
#warning "Compiling with nvcc.\n"
#define HOST_DEVICE __host__ __device__
//#define HOST __host__
#define DEVICE __device__
#define GLOBAL __global__
#include "cuda_runtime.h"
#else
#define HOST_DEVICE
//#define HOST
#define DEVICE
#define GLOBAL
#endif

#include <iostream>
#include <cmath>
#include <cstdint>
#include <cinttypes>
#include <unistd.h>

#define ABS(num) ((num) >= 0 ? (num) : -(num))

#define MILLION 1'000'000
#define BILLION 1'000'000'000
#define TRILLION 1'000'000'000'000
#define QUADRILLION 1'000'000'000'000'000

namespace diff {
    template <typename T>
    concept numeric = requires (T a, T b) {
        {a + b} -> std::same_as<T>;
        {a + b} -> std::same_as<T>;
        {a * b} -> std::same_as<T>;
        {a / b} -> std::same_as<T>;
    };
    template <typename C>
    concept callable = requires (C a) {
        {a()};
    };
    template <typename C, typename T>
    concept callret = requires (C a, T t) {
        {a(t)} -> std::same_as<T>;
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
        [[nodiscard]] char *binrepr(const void *data, uint64_t size) {
            if (!data || !size)
                return nullptr;
            char *str = new char[size*8];
            char *sptr = str;
            const char *ptr = (char*) data;
            while (size --> 0) {
                *sptr++ = ((128 & *ptr) >> 7) + 48;
                *sptr++ = ((64 & *ptr) >> 6) + 48;
                *sptr++ = ((32 & *ptr) >> 5) + 48;
                *sptr++ = ((16 & *ptr) >> 4) + 48;
                *sptr++ = ((8 & *ptr) >> 3) + 48;
                *sptr++ = ((4 & *ptr) >> 2) + 48;
                *sptr++ = ((2 & *ptr) >> 1) + 48;
                *sptr++ = (1 & *ptr++) + 48;
            }
            return str;
        }
    }
    template <numeric T>
    HOST_DEVICE T abs(const T &num) {
        return num >= 0 ? num : -num;
    }
    template <numeric T, callret<T> F>
    HOST_DEVICE T integrate_trap(const F &f, const T &a, const T &b, uint64_t num = 1'000'000) {
        T dx = (b - a)/num;
        T dx_o2 = dx/2;
        uint64_t j = 1;
        T f_i = f(a);
        T f_j;
        T res = 0;
        while (j <= num) {
            res += dx_o2*(f_i + (f_j = f(a + j++*dx)));
            f_i = f_j;
        }
        return res;
    }
    template <numeric T, callret<T> F>
    HOST_DEVICE T trapquad_recurse(const F &f, const T &a, const T &b, const T &tol, const T &f_a, const T &f_b) {
        T f_aPf_b = f_a + f_b;
        T estimate = ((b - a)/2)*f_aPf_b;
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_aPf_b + 2*f_m);
        if (abs(res_or_tol - estimate) > tol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T, callret<T> F>
    HOST_DEVICE T integrate_trapquad(const F &f, const T &a, const T &b, const T &tol = 0.0625l/1024.0l) {
        T f_a = f(a);
        T f_b = f(b);
        T f_aPf_b = f_a + f_b;
        T estimate = ((b - a)/2)*f_aPf_b;
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_aPf_b + 2*f_m);
        if (abs(res_or_tol - estimate) > tol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T, callret<T> F, uint64_t max_depth = 256>
    HOST_DEVICE T trapquad(const F &f, T a, T b, T tol = 0.0625l/1024.0l, uint64_t *mdepth = nullptr) {
        if (mdepth)
            *mdepth = 0;
        uint64_t tree[(max_depth % 64) ? max_depth/64 + 1 : max_depth/64]{}; // fixed at compile-time
        const T global_a = a;
        const T global_b = b;
        T m = (a + b)/2;
        T f_a;// = f(a);
        T f_m;// = f(m);
        T f_b;// = f(b);
        T res = 0;
        T ires; // intermediate result
        T estimate;
        T error;
        T dx = b - a;
        uint64_t depth = 0;
        bool left = true;
        T numdx;
        T whole;
        char *s;
        do {
            start:
            f_a = f(a);
            f_m = f(m);
            f_b = f(b);
            estimate = (dx/2)*(f_a + f_b);
            ires = (dx/4)*(f_a + 2*f_m + f_b);
            if constexpr (max_depth)
                if (depth == max_depth)
                    goto other;
            if (abs(estimate - ires) > tol) {
                tol /= 2;
                dx /= 2;
                left = true;
                b = m;
                m = (a + b)/2; // simplify later
                ++depth;
                // goto start;
            } else {
                other:
                res += ires;
                if (left) {
                    a += dx;
                    m += dx;
                    b += dx;
                    left = false;
                    if constexpr (max_depth)
                        tree[depth/64] |= (1 << (depth % 64));
                    // goto start;
                } else {
                    if (mdepth)
                        if (depth > *mdepth)
                            *mdepth = depth;
                    if (b + dx > global_b)
                        return res;
                    if constexpr (max_depth) {
                        do {
                            tree[depth/64] &= ~(1 << (depth % 64));
                            if (!--depth)
                                return res;
                            a -= dx;
                            dx *= 2;
                            tol *= 2;
                        } while (tree[depth/64] & (1 << (depth % 64)));
                        tree[depth/64] |= (1 << (depth % 64));
                    } else {
                        do {
                            if (!--depth)
                                return res;
                            a -= dx;
                            dx *= 2;
                            tol *= 2;
                            numdx = (b - global_a)/dx;
                            // printf("numdx: %lf\n", numdx);
                            whole = std::round(numdx);
                        } while (std::abs(whole - numdx) > 0.25);
                    }
                    a += dx;
                    // b += dx/2;
                    b += dx;
                    // if (b > global_b)
                    //     return res;
                    m = (a + b)/2;
                    left = false;
                }
            }
            // printf("Depth: %" PRIu64 ", a: %lf, b: %lf, dx: %lf\n", depth, a, b, dx);
            // sleep(1);
        } while (depth);
        return res;
    }
    template <numeric T, callret<T> F>
    HOST_DEVICE T integrate_simpson(const F &f, T a, T b, uint64_t num) {
        T dx = (b - a)/num;
        uint64_t i = 1;
        T offset;
        T a_p_dxo2 = a + dx/2;
        T res = 0.5*f(a) + 2*f(a_p_dxo2);
        while (i < num) {
            offset = i++*dx;
            res += f(a + offset) + 2*f(a_p_dxo2 + offset);
        }
        res += 0.5*f(b);
        res *= dx/3;
        return res;
    }
}
#endif
