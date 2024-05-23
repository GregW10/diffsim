#ifndef DIFFBASE_HPP
#define DIFFBASE_HPP

#if !defined(__linux__) && !defined(__APPLE__)
#error OS not supported
#endif

#include "../glib/misc/gregstack.hpp"
#include <climits>

#if (CHAR_BIT != 8)
#error A byte has to have 8 bits.
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

#ifdef __CUDACC__
#define PI 3.14159265358979323846264338327950288419716939937510582097494459
#else
#define PI 3.14159265358979323846264338327950288419716939937510582097494459l
#endif

#define EMPTY

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
    HOST_DEVICE T trapquad_recurse(const F &f, const T &a, const T &b, const T &tol, const T &ptol,
                                   const T &f_a, const T &f_b) {
        T f_aPf_b = f_a + f_b;
        T estimate = ((b - a)/2)*f_aPf_b;
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_aPf_b + 2*f_m);
        if (abs(res_or_tol - estimate) > tol ||
            std::abs((a - midpoint)*(f_b - f_m) + (f_m - f_a)*(b - midpoint)) <= ptol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, ptol/2, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, ptol/2, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T, callret<T> F>
    HOST_DEVICE T integrate_trapquad(const F &f, const T &a, const T &b, const T &tol = 0.0625l/1024.0l,
                                     const T &ptol = 0.0625l/4096.0l) {
        T f_a = f(a);
        T f_b = f(b);
        T f_aPf_b = f_a + f_b;
        T estimate = ((b - a)/2)*f_aPf_b;
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_aPf_b + 2*f_m);
        if (abs(res_or_tol - estimate) > tol ||
            std::abs((a - midpoint)*(f_b - f_m) + (f_m - f_a)*(b - midpoint)) <= ptol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, ptol/2, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, ptol/2, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T, callret<T> F, uint64_t max_depth = 256, bool use_ptol = true>
    requires (std::is_floating_point<T>::value)
    HOST_DEVICE T trapquad(const F &f, T a, T b, T tol = 0.0625l/1024.0l, T ptol = 0, uint64_t *mdepth = nullptr) {
        if (tol < 0)
            return std::numeric_limits<T>::quiet_NaN();
        if constexpr (use_ptol)
            if (ptol < 0)
                return std::numeric_limits<T>::quiet_NaN();
        if (mdepth)
            *mdepth = 0;
        struct func_vals {
            T fm;
            T fb;
        } fmfb_val;
        gtd::stack<func_vals> stack;
        uint64_t tree[(max_depth % 64) ? max_depth/64 + 1 : max_depth/64]{}; // fixed at compile-time
        const T global_a = a;
        const T global_b = b;
        T m = (a + b)/2;
        T f_a = f(a);
        T f_m = f(m);
        T f_b = f(b);
        T res = 0;
        T ires; // intermediate result
        T estimate;
        // T error;
        T dx = b - a;
        T dxo2 = dx/2;
        T dxo4 = dx/4;
        uint64_t depth = 0;
        bool left = false;
        T numdx;
        T whole;
        T f_aPf_b;
        do {
            start:
            // f_a = f(a);
            // f_m = f(m);
            // f_b = f(b);
            f_aPf_b = f_a + f_b;
            estimate = dxo2*f_aPf_b;
            ires = dxo4*(f_aPf_b + 2*f_m);
            if constexpr (max_depth)
                if (depth == max_depth)
                    goto other;
            if constexpr (use_ptol)
                if (std::abs(f_m - f_a) <= ptol && std::abs(f_b - f_m) <= ptol) // < 2*ptol
                    goto divide;
            if (abs(estimate - ires) > tol) {
                divide:
                printf("Split into two.\n");
                dx = dxo2;
                dxo2 = dxo4;
                dxo4 /= 2;
                // dx /= 2;
                tol /= 2;
                if constexpr (use_ptol)
                    ptol /= 2;
                left = true;
                b = m;
                // m = (a + b)/2; // simplify later
                m -= dxo2;
                fmfb_val.fm = f_m;
                fmfb_val.fb = f_b;
                stack.push(fmfb_val);
                f_b = f_m;
                f_m = f(m);
                ++depth;
                // goto start;
            } else {
                other:
                // printf("a: %lf, b: %lf\n", a, b);
                res += ires;
                if (left) {
                    printf("else -> if (left)\n");
                    a = b;
                    // a += dx;
                    m += dx;
                    b += dx;
                    f_a = f_b;
                    printf("Size of stack: %" PRIu64 "\n", stack.size());
                    f_m = f(m);
                    f_b = f(b);
                    // printf("f(m) = %lf, top: %lf\n", f_m, stack.top().fm);
                    printf("f(b) = %lf, top: %lf\n", f_b, stack.top().fb);
                    fmfb_val = stack.pop_r();
                    f_m = fmfb_val.fm;
                    f_b = fmfb_val.fb;
                    // f_b = stack.top().fb;
                    // stack.top().fm = f_m;
                    left = false;
                    if constexpr (max_depth)
                        tree[depth/64] |= (1 << (depth % 64));
                    // goto start;
                } else {
                    printf("else -> else\n");
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
                            dxo4 = dxo2;
                            dxo2 = dx;
                            dx *= 2;
                            tol *= 2;
                            if constexpr (use_ptol)
                                ptol *= 2;
                            f_m = stack.pop_r().fb;
                            printf("Stack size: %" PRIu64 ", depth: %" PRIu64 "\n", stack.size(), depth);
                        } while (tree[depth/64] & (1 << (depth % 64)));
                        tree[depth/64] |= (1 << (depth % 64));
                    } else {
                        do {
                            if (!--depth)
                                return res;
                            // a -= dx;
                            dxo4 = dxo2;
                            dxo2 = dx;
                            dx *= 2;
                            tol *= 2;
                            if constexpr (use_ptol)
                                ptol *= 2;
                            f_m = stack.pop_r().fb;
                            numdx = (b - global_a)/dx;
                            // printf("numdx: %lf\n", numdx);
                            whole = std::round(numdx);
                        } while (std::abs(whole - numdx) > 0.25);
                    }
                    // a += dx;
                    a = b;
                    // b += dx/2;
                    b += dx;
                    // if (b > global_b)
                    //     return res;
                    m = (a + b)/2;
                    f_a = f_b;
                    f_b = f_m;
                    f_m = f(m);
                    fmfb_val.fm = f_m;
                    fmfb_val.fb = f_b;
                    stack.push(fmfb_val);
                    // f_b = f(b);
                    left = false;
                }
            }
            // printf("Depth: %" PRIu64 ", a: %lf, b: %lf, dx: %lf\n", depth, a, b, dx);
            // sleep(1);
        } while (depth);
        std::cout << "Size of stack: " << stack.size() << std::endl;
        return res;
    }
#define NO_MAX_DEPTH ((uint64_t) -1)
    template <numeric T, callret<T> F>
    requires (std::is_floating_point<T>::value)
    HOST_DEVICE T trapquad_n(const F &f, T a, T b, T tol = 0.0625l/1024.0l, T ptol = 1/(1024.0l*1024.0l),
        uint64_t *mdepth = nullptr) {
        if (b <= a || tol <= 0)
            return std::numeric_limits<T>::quiet_NaN();
        gtd::stack<T> stack;
        T dx = b - a;
        T dxo2 = dx/2;
        T dxo4 = dx/4;
        T m = (a + b)/2;
        T fa = f(a);
        T fm = f(m);
        T fb = f(b);
        T faPfb;
        T estimate;
        T ires;
        T res = 0;
        bool left = false;
#define CROSS std::abs(dxo2*((fm - fb) + (fm - fa)))
#define TRAPQUAD_LOOP(par0, par00, par1, par2, par3, lab1, lab2) \
        while (1) { \
            par1 \
            faPfb = fa + fb; \
            estimate = dxo2*faPfb; \
            ires = dxo4*(faPfb + 2*fm); \
            par2 \
            if (std::abs(ires - estimate) > tol) { \
                lab1: \
                dx = dxo2; \
                dxo2 = dxo4; \
                dxo4 /= 2; \
                tol /= 2; \
                par0 \
                b = m; \
                m -= dxo2; \
                stack.push(fb); \
                fb = fm; \
                fm = f(m); \
                left = true; \
            } \
            else { \
                lab2: \
                res += ires; \
                if (left) { \
                    a = b; \
                    m += dx; \
                    b += dx; \
                    fa = fb; \
                    fm = f(m); \
                    fb = stack.top(); \
                    left = false; \
                } \
                else { \
                    if (!stack) \
                        return res; \
                    par3 \
                    while (stack.top() == fb) { \
                        m = a; \
                        a -= dx; \
                        dxo4 = dxo2; \
                        dxo2 = dx; \
                        dx *= 2; \
                        tol *= 2; \
                        par00 \
                        stack.pop(); \
                        if (!stack) \
                            return res; \
                    } \
                    a = b; \
                    m += dx; \
                    b += dx; \
                    fa = fb; \
                    fm = f(m); \
                    fb = stack.top(); \
                    left = false; \
                } \
            } \
        }
        if (ptol >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;, EMPTY,
                                  if (CROSS <= ptol) {goto divide_1;},
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_1, other_1)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                                  if ((_s = stack.size()) == max_depth) {ires = dxo4*(faPfb + 2*fm); goto other_2;},
                                  if (CROSS <= ptol) {goto divide_2;},
                                  if (_s > *mdepth) {*mdepth = _s;}, divide_2, other_2)
                }
            } else {
                TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                              EMPTY, if (CROSS <= ptol) {goto divide_3;},
                              EMPTY, divide_3, other_3)
            }
        } else {
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    TRAPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY,
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_4, other_4)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    TRAPQUAD_LOOP(EMPTY, EMPTY,
                                  if ((_s = stack.size()) == max_depth) {ires = dxo4*(faPfb + 2*fm); goto other_5;},
                                  EMPTY, if (_s > *mdepth) {*mdepth = _s;}, divide_5, other_5)
                }
            } else {
                TRAPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_6, other_6)
            }
        }
#undef TRAPQUAD_LOOP
#undef CROSS
        return res; // for completeness, but would never be reached
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
