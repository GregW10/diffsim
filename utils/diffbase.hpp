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
// #include <fstream>

#define ABS(num) ((num) >= 0 ? (num) : -(num))

#define MILLION 1'000'000
#define BILLION 1'000'000'000
#define TRILLION 1'000'000'000'000
#define QUADRILLION 1'000'000'000'000'000

#ifdef __CUDACC__
#define PI 3.14159265358979323846264338327950288419716939937510582097494459
#define LIGHT_SPEED 299'792'458.0 // m/s
#define PERMITTIVITY 0.0000000000088541878188 // F/m
#define PERMEABILITY 0.00000125663706127 // F/m
#else
#define PI 3.14159265358979323846264338327950288419716939937510582097494459l
#define LIGHT_SPEED 299'792'458.0l // m/s
#define PERMITTIVITY 0.0000000000088541878188l // F/m
#define PERMEABILITY 0.00000125663706127 // F/m
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
    template <typename C, typename T>
    concept calldblret = requires (C a, T t1, T t2) {
        {a(t1, t2)} -> std::same_as<T>;
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
        template <numeric T>
        HOST_DEVICE T abs(const T &num) {
            return num >= 0 ? num : -num;
        }
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
#define NO_MAX_DEPTH ((uint64_t) -1)
    template <numeric T, callret<T> F, bool prog = false>
    requires (std::is_floating_point<T>::value)
    HOST_DEVICE T trapquad(const F &f, T a, T b, T tol = 0.0625l/1024.0l, T ptol = 1/(1024.0l*1024.0l),
        uint64_t *mdepth = nullptr) {
        if (b < a || tol <= 0)
            return std::numeric_limits<T>::quiet_NaN();
        if (a == b) { // zero width = zero area, so might as well save some computation below
            if (mdepth)
                *mdepth = 0;
            return 0;
        }
        gtd::stack<T> stack;
        T global_range = b - a;
        T global_a = a;
        T dx = global_range;
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
        uint64_t bitmask;
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
                if ((stack.size()) > tree.size()*64) { \
                    tree.push(); \
                } \
            } \
            else { \
                lab2: \
                res += ires; \
                if constexpr (prog) { \
                    printf("%.5Lf%% complete\r", ((((long double) b) - global_a)/global_range)*100); \
                    fflush(stdout); \
                } \
                if (left) { \
                    a = b; \
                    m += dx; \
                    b += dx; \
                    fa = fb; \
                    fm = f(m); \
                    fb = stack.top(); \
                    left = false; \
                    tree.top() |= (1ull << ((stack.size() - 1) % 64)); \
                } \
                else { \
                    if (!stack) { \
                        if constexpr (prog) \
                            putchar('\n'); \
                        return res; \
                    } \
                    par3 \
                    do { \
                        stack.pop(); \
                        if (!stack) { \
                            if constexpr (prog) \
                                putchar('\n'); \
                            return res; \
                        } \
                        tree.top() &= ~(1ull << (stack.size() % 64)); \
                        m = a; \
                        a -= dx; \
                        dxo4 = dxo2; \
                        dxo2 = dx; \
                        dx *= 2; \
                        tol *= 2; \
                        par00 \
                        bitmask = (1ull << ((stack.size() - 1) % 64)); \
                        if (bitmask == 9'223'372'036'854'775'808u) \
                            tree.pop(); \
                    } while ((tree.top() & bitmask) > 0); \
                    a = b; \
                    m += dx; \
                    b += dx; \
                    fa = fb; \
                    fm = f(m); \
                    fb = stack.top(); \
                    left = false; \
                    tree.top() |= bitmask; \
                } \
            } \
        }
        if (ptol >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                    TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;, EMPTY,
                                  if (CROSS <= ptol) {goto divide_1;},
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_1, other_1)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                                  if ((_s = stack.size()) == max_depth) {ires = dxo4*(faPfb + 2*fm); goto other_2;},
                                  if (CROSS <= ptol) {goto divide_2;},
                                  if (_s > *mdepth) {*mdepth = _s;}, divide_2, other_2)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                TRAPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                              EMPTY, if (CROSS <= ptol) {goto divide_3;},
                              EMPTY, divide_3, other_3)
            }
        } else {
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8};
                    TRAPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY,
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_4, other_4)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    TRAPQUAD_LOOP(EMPTY, EMPTY,
                                  if ((_s = stack.size()) == max_depth) {ires = dxo4*(faPfb + 2*fm); goto other_5;},
                                  EMPTY, if (_s > *mdepth) {*mdepth = _s;}, divide_5, other_5)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                TRAPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_6, other_6)
            }
        }
#undef TRAPQUAD_LOOP
#undef CROSS
        if constexpr (prog)
            putchar('\n');
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
    template <numeric T, callret<T> F, bool prog = false>
    requires (std::is_floating_point<T>::value)
    HOST_DEVICE T simpquad(const F &f, T a, T b, T abstol = 0.0625l/1024.0l, T reltol = 0.0625l/1024.0l,
        T ptol = 1/(1024.0l*1024.0l), uint64_t *mdepth = nullptr) {
        if (b < a || abstol <= 0 || reltol < 0 || reltol >= 1)
            return std::numeric_limits<T>::quiet_NaN();
        if (a == b) {// zero width = zero area, so might as well save some computation below
            if (mdepth)
                *mdepth = 0;
            return 0;
        }
        struct fm3fb_val {
            T fm3;
            T fb;
        } fm3fb;
        gtd::stack<fm3fb_val> stack;
        T global_range = b - a;
        T global_a = a;
        T dx = global_range;
        T dxo2 = dx/2;
        T dxo4 = dx/4;
        T dxo6 = dx/6;
        T dxo12 = dx/12;
        T m2 = (a + b)/2;
        T m1 = m2 - dxo4;
        T m3 = m2 + dxo4;
        T fa = f(a);
        T fm1 = f(m1);
        T fm2 = f(m2);
        T fm3 = f(m3);
        T fb = f(b);
        T faPfb;
        T estimate;
        T ires;
        T res = 0;
        bool left = false;
        uint64_t bitmask;
        T absval;
#define CROSS std::abs(dxo2*((fm2 - fb) + (fm2 - fa)))
#define SIMPQUAD_LOOP(par0, par00, par1, par2, par3, lab1, lab2) \
        while (1) { \
            par1 \
            faPfb = fa + fb; \
            estimate = dxo6*(faPfb + 4*fm2); \
            ires = dxo12*(faPfb + 4*fm1 + 2*fm2 + 4*fm3); \
            par2 \
            if ((absval = std::abs(ires - estimate)) > abstol && absval > std::abs(reltol*ires)) { \
                lab1: \
                dx = dxo2; \
                dxo2 = dxo4; \
                dxo4 /= 2; \
                dxo6 = dxo12; \
                dxo12 /= 2; \
                abstol /= 2; \
                par0 \
                b = m2; \
                m2 -= dxo2; \
                m1 -= dxo4; \
                m3 = m2 + dxo4; \
                fm3fb.fm3 = fm3; \
                fm3fb.fb = fb; \
                stack.push(fm3fb); \
                fb = fm2; \
                fm2 = fm1; \
                fm1 = f(m1); \
                fm3 = f(m3); \
                left = true; \
                if ((stack.size()) > tree.size()*64) { \
                    tree.push(); \
                } \
            } \
            else { \
                lab2: \
                res += ires; \
                if constexpr (prog) { \
                    printf("%.5Lf%% complete\r", ((((long double) b) - global_a)/global_range)*100); \
                    fflush(stdout); \
                } \
                if (left) { \
                    a = b; \
                    m1 += dx; \
                    m2 += dx; \
                    m3 += dx; \
                    b += dx; \
                    fa = fb; \
                    fm1 = f(m1); \
                    fm3 = f(m3); \
                    fm3fb = stack.top(); \
                    fm2 = fm3fb.fm3; \
                    fb = fm3fb.fb; \
                    left = false; \
                    tree.top() |= (1ull << ((stack.size() - 1) % 64)); \
                } \
                else { \
                    if (!stack) { \
                        if constexpr (prog) \
                            putchar('\n'); \
                        return res; \
                    } \
                    par3 \
                    do { \
                        stack.pop(); \
                        if (!stack) { \
                            if constexpr (prog) \
                                putchar('\n'); \
                            return res; \
                        } \
                        tree.top() &= ~(1ull << (stack.size() % 64)); \
                        m2 = a; \
                        a -= dx; \
                        m1 = a + dxo2; \
                        m3 = m1 + dx; \
                        dxo4 = dxo2; \
                        dxo2 = dx; \
                        dx *= 2; \
                        dxo12 = dxo6; \
                        dxo6 *= 2; \
                        abstol *= 2; \
                        par00 \
                        bitmask = (1ull << ((stack.size() - 1) % 64)); \
                        if (bitmask == 9'223'372'036'854'775'808u) \
                            tree.pop(); \
                    } while ((tree.top() & bitmask) > 0); \
                    a = b; \
                    m1 += dx; \
                    m2 += dx; \
                    m3 += dx; \
                    b += dx; \
                    fa = fb; \
                    fm1 = f(m1); \
                    fm3 = f(m3); \
                    fm3fb = stack.top(); \
                    fm2 = fm3fb.fm3; \
                    fb = fm3fb.fb; \
                    left = false; \
                    tree.top() |= bitmask; \
                } \
            } \
        }
        if (ptol >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                    SIMPQUAD_LOOP(ptol /= 2;, ptol *= 2;, EMPTY,
                                  if (CROSS <= ptol) {goto divide_1;},
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_1, other_1)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    SIMPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                                  if ((_s = stack.size()) == max_depth)
                                  {ires = dxo12*(fa + 4*fm1 + 2*fm2 + 4*fm3 + fb); goto other_2;},
                                  if (CROSS <= ptol) {goto divide_2;},
                                  if (_s > *mdepth) {*mdepth = _s;}, divide_2, other_2)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                SIMPQUAD_LOOP(ptol /= 2;, ptol *= 2;,
                              EMPTY, if (CROSS <= ptol) {goto divide_3;},
                              EMPTY, divide_3, other_3)
            }
        } else {
            if (mdepth) {
                uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8};
                    SIMPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY,
                                  if ((_s = stack.size()) > *mdepth) {*mdepth = _s;}, divide_4, other_4)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    SIMPQUAD_LOOP(EMPTY, EMPTY,
                                  if ((_s = stack.size()) == max_depth)
                                  {ires = dxo12*(fa + 4*fm1 + 2*fm2 + 4*fm3 + fb); goto other_5;},
                                  EMPTY, if (_s > *mdepth) {*mdepth = _s;}, divide_5, other_5)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                SIMPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_6, other_6)
            }
        }
#undef SIMPQUAD_LOOP
#undef CROSS
        if constexpr (prog)
            putchar('\n');
        return res; // for completeness, but would never be reached
    }
    template <numeric T, calldblret<T> F, callret<T> GH, bool prog = false>
    requires (std::is_floating_point<T>::value)
    HOST_DEVICE T simpdblquad(const F &f,
                              T ya,
                              T yb,
                              const GH &gfunc,
                              const GH &hfunc,
                              T abstol_y = 0.0625l/1024.0l,
                              T reltol_y = 0.0625l/1024.0l,
                              T ptol_y = 1/(1024.0l*1024.0l),
                              uint64_t *mdepth_y = nullptr,
                              T abstol_x = 0.0625l/1024.0l,
                              T reltol_x = 0.0625l/1024.0l,
                              T ptol_x = 1/(1024.0l*1024.0l),
                              uint64_t *mdepth_x = nullptr) {
        if (yb <= ya || abstol_y <= 0)
            return std::numeric_limits<T>::quiet_NaN();
        if (ya == yb) // zero width = zero area, so might as well save some computation below
            return 0;
        struct fym3fyb_val {
            T fym3;
            T fyb;
        } fym3fyb;
#define XDEPTH if (maxx_depth > *mdepth_x) *mdepth_x = maxx_depth; maxx_depth = org_mxd;
#define CROSS std::abs(dyo2*((fym2 - fyb) + (fym2 - fya)))
        gtd::stack<fym3fyb_val> stack;
        T global_yrange = yb - ya;
        T global_ya = ya;
        T dy = global_yrange;
        T dyo2 = dy/2;
        T dyo4 = dy/4;
        T dyo6 = dy/6;
        T dyo12 = dy/12;
        T ym2 = (ya + yb)/2;
        T ym1 = ym2 - dyo4;
        T ym3 = ym2 + dyo4;
        T y;
        auto fy_lam = [&f, &y](T x){return f(x, y);};
        T fyaPfyb;
        T estimate;
        T ires;
        T res = 0;
        bool left = false;
        uint64_t bitmask;
        T absval;
        T fya;
        T fym1;
        T fym2;
        T fym3;
        T fyb;
        uint64_t *mptr;
#define SIMPDBLQUAD_MOST(par0, par00, par1, par2, par3, lab1, lab2, after) \
        y = ya; \
        fya = simpquad(fy_lam, gfunc(ya), hfunc(ya), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym1; \
        fym1 = simpquad(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym2; \
        fym2 = simpquad(fy_lam, gfunc(ym2), hfunc(ym2), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym3; \
        fym3 = simpquad(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr); after \
        y = yb; \
        fyb = simpquad(fy_lam, gfunc(yb), hfunc(yb), abstol_x, reltol_x, ptol_x, mptr); after \
        while (1) { \
            par1 \
            fyaPfyb = fya + fyb; \
            estimate = dyo6*(fyaPfyb + 4*fym2); \
            ires = dyo12*(fyaPfyb + 4*fym1 + 2*fym2 + 4*fym3); \
            par2 \
            if ((absval = std::abs(ires - estimate)) > abstol_y && absval > reltol_y*ires) { \
                lab1: \
                dy = dyo2; \
                dyo2 = dyo4; \
                dyo4 /= 2; \
                dyo6 = dyo12; \
                dyo12 /= 2; \
                abstol_y /= 2; \
                par0 \
                yb = ym2; \
                ym2 -= dyo2; \
                ym1 -= dyo4; \
                ym3 = ym2 + dyo4; \
                fym3fyb.fym3 = fym3; \
                fym3fyb.fyb = fyb; \
                stack.push(fym3fyb); \
                fyb = fym2; \
                fym2 = fym1; \
                y = ym1; \
                fym1 = simpquad(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr); after \
                y = ym3; \
                fym3 = simpquad(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr); after \
                left = true; \
                if ((stack.size()) > tree.size()*64) { \
                    tree.push(); \
                } \
            } \
            else { \
                lab2: \
                if constexpr (prog) { \
                    printf("%.5Lf%% complete\r", ((((long double) yb) - global_ya)/global_yrange)*100); \
                    fflush(stdout); \
                } \
                res += ires; \
                if (left) { \
                    ya = yb; \
                    ym1 += dy; \
                    ym2 += dy; \
                    ym3 += dy; \
                    yb += dy; \
                    fya = fyb; \
                    y = ym1; \
                    fym1 = simpquad(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr); after \
                    y = ym3; \
                    fym3 = simpquad(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr); after \
                    fym3fyb = stack.top(); \
                    fym2 = fym3fyb.fym3; \
                    fyb = fym3fyb.fyb; \
                    left = false; \
                    tree.top() |= (1ull << ((stack.size() - 1) % 64)); \
                } \
                else { \
                    if (!stack) { \
                        if constexpr (prog) \
                            putchar('\n'); \
                        return res; \
                    } \
                    par3 \
                    do { \
                        stack.pop(); \
                        if (!stack) { \
                            if constexpr (prog) \
                                putchar('\n'); \
                            return res; \
                        } \
                        tree.top() &= ~(1ull << (stack.size() % 64)); \
                        ym2 = ya; \
                        ya -= dy; \
                        ym1 = ya + dyo2; \
                        ym3 = ym1 + dy; \
                        dyo4 = dyo2; \
                        dyo2 = dy; \
                        dy *= 2; \
                        dyo12 = dyo6; \
                        dyo6 *= 2; \
                        abstol_y *= 2; \
                        par00 \
                        bitmask = (1ull << ((stack.size() - 1) % 64)); \
                        if (bitmask == 9'223'372'036'854'775'808u) \
                            tree.pop(); \
                        /*printf("stack.size(): %" PRIu64 ", tree.top(): %" PRIu64 ", bitmask: %" PRIu64 "\n",*/ \
                        /*stack.size(), tree.top(), bitmask); */\
                    } while ((tree.top() & bitmask) > 0); \
                    ya = yb; \
                    ym1 += dy; \
                    ym2 += dy; \
                    ym3 += dy; \
                    yb += dy; \
                    fya = fyb; \
                    y = ym1; \
                    fym1 = simpquad(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr); after \
                    y = ym3; \
                    fym3 = simpquad(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr); after \
                    fym3fyb = stack.top(); \
                    fym2 = fym3fyb.fym3; \
                    fyb = fym3fyb.fyb; \
                    left = false; \
                    tree.top() |= bitmask; \
                } \
            } \
        }
        if (ptol_y >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            if (mdepth_y) {
                uint64_t _s;
                if (*mdepth_y == NO_MAX_DEPTH) {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                        SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;, EMPTY,
                                         if (CROSS <= ptol_y) {goto divide_11;},
                                         if ((_s = stack.size()) > *mdepth_y) {*mdepth_y = _s;},
                                         divide_11, other_11, XDEPTH)
                    } else {
                        mptr = nullptr;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                        SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;, EMPTY,
                                         if (CROSS <= ptol_y) {goto divide_12;},
                                         if ((_s = stack.size()) > *mdepth_y) {*mdepth_y = _s;},
                                         divide_12, other_12, EMPTY)
                    }
                } else {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;,
                                         if ((_s = stack.size()) == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_21;},
                                         if (CROSS <= ptol_y) {goto divide_21;},
                                         if (_s > *mdepth_y) {*mdepth_y = _s;}, divide_21, other_21, XDEPTH)
                    } else {
                        mptr = nullptr;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;,
                                         if ((_s = stack.size()) == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_22;},
                                         if (CROSS <= ptol_y) {goto divide_22;},
                                         if (_s > *mdepth_y) {*mdepth_y = _s;}, divide_22, other_22, EMPTY)
                    }
                }
            } else {
                if (mdepth_x) {
                    uint64_t maxx_depth = *mdepth_x;
                    const uint64_t org_mxd = *mdepth_x;
                    *mdepth_x = 0;
                    mptr = &maxx_depth;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;,
                                     EMPTY, if (CROSS <= ptol_y) {goto divide_31;},
                                     EMPTY, divide_31, other_31, XDEPTH)
                } else {
                    mptr = nullptr;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(ptol_y /= 2;, ptol_y *= 2;,
                                     EMPTY, if (CROSS <= ptol_y) {goto divide_32;},
                                     EMPTY, divide_32, other_32, EMPTY)
                }
            }
        } else {
            if (mdepth_y) {
                uint64_t _s;
                if (*mdepth_y == NO_MAX_DEPTH) {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY,
                                         if ((_s = stack.size()) > *mdepth_y) {*mdepth_y = _s;},
                                         divide_41, other_41, XDEPTH)
                    } else {
                        mptr = nullptr;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY,
                                         if ((_s = stack.size()) > *mdepth_y) {*mdepth_y = _s;},
                                         divide_42, other_42, EMPTY)
                    }
                } else {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY,
                                         if ((_s = stack.size()) == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_51;},
                                         EMPTY, if (_s > *mdepth_y) {*mdepth_y = _s;}, divide_51, other_51, XDEPTH)
                    } else {
                        mptr = nullptr;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY,
                                         if ((_s = stack.size()) == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_52;},
                                         EMPTY, if (_s > *mdepth_y) {*mdepth_y = _s;}, divide_52, other_52, EMPTY)
                    }
                }
            } else {
                if (mdepth_x) {
                    uint64_t maxx_depth = *mdepth_x;
                    const uint64_t org_mxd = *mdepth_x;
                    *mdepth_x = 0;
                    mptr = &maxx_depth;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_61, other_61, XDEPTH)
                } else {
                    mptr = nullptr;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_62, other_62, EMPTY)
                }
            }
        }
#undef SIMPDBLQUAD_MOST
#undef CROSS
        if constexpr (prog)
            putchar('\n');
        return res; // for completeness, but would never be reached
    }
}
#endif
