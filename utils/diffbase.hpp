#ifndef DIFFBASE_HPP
#define DIFFBASE_HPP

#if !defined(__linux__) && !defined(__APPLE__)
#error OS not supported
#endif

#include "../glib/misc/gregcomplex.hpp"
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

#define EMPTY

namespace diff {
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
        template <gtd::numeric T>
        HOST_DEVICE T abs(const T &num) {
            return num >= 0 ? num : -num;
        }
    }
    template <gtd::numeric T>
    HOST_DEVICE T Efield_to_intensity(const T &Efield) {
        return LIGHT_SPEED*PERMITTIVITY*Efield*Efield;
    }
    template <gtd::numeric T>
    HOST_DEVICE T Efield_to_intensity(const gtd::complex<T> &Efield) {
        return LIGHT_SPEED*PERMITTIVITY*static_cast<T>(Efield*Efield.conj()); // change to Efield.mag_sq()
    }
    template <gtd::numeric T> // impossible to recover phase information, so must return scalar (real)
    HOST_DEVICE T intensity_to_Efield(const T &intensity) {
        return std::sqrt(intensity/(LIGHT_SPEED*PERMITTIVITY));
    }
    template <gtd::numeric T>
    HOST_DEVICE T E0_to_intensity(const T &E0) {
        return 0.5l*LIGHT_SPEED*PERMITTIVITY*E0*E0;
    }
    template <gtd::numeric T>
    HOST_DEVICE T E0_to_intensity(const gtd::complex<T> &E0) {
        return 0.5l*LIGHT_SPEED*PERMITTIVITY*static_cast<T>(E0*E0.conj()); // change to Efield.mag_sq()
    }
    template <gtd::numeric T>
    HOST_DEVICE T intensity_to_E0(const T &intensity) {
        return std::sqrt((2*intensity)/(LIGHT_SPEED*PERMITTIVITY));
    }
    template <gtd::numeric T, gtd::callret<T> F>
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
    template <gtd::numeric T, gtd::callret<T> F>
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
    template <gtd::numeric T, gtd::callret<T> F>
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
    template <gtd::numeric T, gtd::callret<T> F, bool prog = false>
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
    template <gtd::numeric T, gtd::callret<T> F>
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
    template <gtd::numeric T, gtd::numeric R, gtd::callret<T, R> F, bool prog>
    HOST_DEVICE R simpquad(const F &f,
                           T a,
                           T b,
                           T abstol,
                           T reltol,
                           T ptol,
                           uint64_t *mdepth) {
        if (b < a || abstol <= 0 || reltol < 0 || reltol >= 1)
            throw std::invalid_argument{"Error: upper x-bound must be above lower x-bound, x-tolerance must be "
                                        "positive and x-relative-tolerance must be in [0,1).\n"};
        if (a == b) {// zero width = zero area, so might as well save some computation below
            if (mdepth)
                *mdepth = 0;
            return 0;
        }
        struct fm3fb_val {
            R fm3;
            R fb;
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
        R fa = f(a);
        R fm1 = f(m1);
        R fm2 = f(m2);
        R fm3 = f(m3);
        R fb = f(b);
        R fm1mfm2;
        R fm2mfm3;
        R faPfb;
        R estimate;
        R ires;
        R res = 0;
        bool left = false;
        uint64_t bitmask;
        T absval;
// #define CROSS std::abs(dxo2*((fm2 - fb) + (fm2 - fa)))
#define BT fm1mfm2 = fm1 - fm2; fm2mfm3 = fm2 - fm3;
#define CROSS1 std::abs(dxo4*(fm1mfm2 + (fm1 - fa))) // see if I can save on computation here
#define CROSS2 std::abs(dxo4*(fm2mfm3 - fm1mfm2))
#define CROSS3 std::abs(dxo4*((fm3 - fb) - fm2mfm3))
#define CCROSSr1 std::abs(dxo4*(fm1mfm2._real + (fm1._real - fa._real))) // see if I can save on computation here
#define CCROSSr2 std::abs(dxo4*(fm2mfm3._real - fm1mfm2._real))
#define CCROSSr3 std::abs(dxo4*((fm3._real - fb._real) - fm2mfm3._real))
#define CCROSSi1 std::abs(dxo4*(fm1mfm2._imag + (fm1._imag - fa._imag))) // see if I can save on computation here
#define CCROSSi2 std::abs(dxo4*(fm2mfm3._imag - fm1mfm2._imag))
#define CCROSSi3 std::abs(dxo4*((fm3._imag - fb._imag) - fm2mfm3._imag))
#define CROSS(num) \
        if constexpr (std::same_as<R, gtd::complex<T>>) { \
            if (CCROSSr1 <= ptol && CCROSSr2 <= ptol && CCROSSr3 <= ptol && \
                CCROSSi1 <= ptol && CCROSSi2 <= ptol && CCROSSi3 <= ptol) {goto divide_##num;} \
        } else { /* printf("c1: %lf, c2: %lf, c3: %lf, ptol: %lf\n", CROSS1, CROSS2, CROSS3, ptol); */ \
            if (CROSS1 <= ptol && CROSS2 <= ptol && CROSS3 <= ptol) {goto divide_##num;} \
        }
#define SIMPQUAD_LOOP(par0, par00, par1, par2, par3, lab1, lab2) \
        while (1) { \
            par1 \
            faPfb = fa + fb; \
            estimate = dxo6*(faPfb + 4*fm2); \
            ires = dxo12*(faPfb + 4*fm1 + 2*fm2 + 4*fm3); \
            par2 \
            if ((absval = gtd::abs(ires - estimate)) > abstol && absval > gtd::abs(reltol*ires)) { \
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
                if (stack._size > tree.size()*64) { \
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
                        /* par00 */ \
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
                    par00 \
                } \
            } \
        }
        if (ptol >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            const T org_ptol = ptol;
            if (mdepth) {
                // uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                    SIMPQUAD_LOOP(ptol /= 4;, ptol = org_ptol/(1 << 2*stack._size);, EMPTY, BT CROSS(1),
                                  if (stack._size > *mdepth) {*mdepth = stack._size;}, divide_1, other_1)
                } else {
                    // printf("in here.\n");
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    SIMPQUAD_LOOP(ptol /= 4;, ptol = org_ptol/(1 << 2*stack._size);,
                                  if (stack._size == max_depth)
                                  {ires = dxo12*(fa + 4*fm1 + 2*fm2 + 4*fm3 + fb); goto other_2;}, BT
                                  CROSS(2),
                                  if (stack._size > *mdepth) {*mdepth = stack._size;}, divide_2, other_2)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                SIMPQUAD_LOOP(ptol /= 4;, ptol = org_ptol/(1 << 2*stack._size);,
                              EMPTY, BT CROSS(3),
                              EMPTY, divide_3, other_3)
            }
        } else {
            if (mdepth) {
                // uint64_t _s;
                if (*mdepth == NO_MAX_DEPTH) {
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{8};
                    SIMPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY,
                                  if (stack._size > *mdepth) {*mdepth = stack._size;}, divide_4, other_4)
                } else {
                    uint64_t max_depth = *mdepth;
                    *mdepth = 0;
                    gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                    SIMPQUAD_LOOP(EMPTY, EMPTY,
                                  if (stack._size == max_depth)
                                  {ires = dxo12*(fa + 4*fm1 + 2*fm2 + 4*fm3 + fb); goto other_5;},
                                  EMPTY, if (stack._size > *mdepth) {*mdepth = stack._size;}, divide_5, other_5)
                }
            } else {
                gtd::stack<uint64_t> tree{8};
                SIMPQUAD_LOOP(EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, divide_6, other_6)
            }
        }
#undef SIMPQUAD_LOOP
#undef CROSS
#undef CCROSSi3
#undef CCROSSi2
#undef CCROSSi1
#undef CCROSSr3
#undef CCROSSr2
#undef CCROSSr1
#undef CROSS1
#undef CROSS2
#undef CROSS3
#undef BT
        if constexpr (prog)
            putchar('\n');
        return res; // for completeness, but would never be reached
    }
    template <gtd::numeric T, gtd::numeric R, gtd::calldblret<T, R> F, gtd::callret<T> GH, bool prog>
    HOST_DEVICE R simpdblquad(const F &f,
                              T ya,
                              T yb,
                              const GH &gfunc,
                              const GH &hfunc,
                              T abstol_y,
                              T reltol_y,
                              T ptol_y,
                              uint64_t *mdepth_y,
                              T abstol_x,
                              T reltol_x,
                              T ptol_x,
                              uint64_t *mdepth_x) {
        if (yb <= ya || abstol_y <= 0 || reltol_y < 0 || reltol_y >= 1)
            throw std::invalid_argument{"Error: upper y-bound must be above lower y-bound, y-tolerance must be "
                                        "positive and y-relative-tolerance must be in [0,1).\n"};
        if (ya == yb) // zero width = zero area, so might as well save some computation below
            return 0;
        struct fym3fyb_val {
            R fym3;
            R fyb;
        } fym3fyb;
#define XDEPTH if (maxx_depth > *mdepth_x) *mdepth_x = maxx_depth; maxx_depth = org_mxd;
// #define CROSS std::abs(dyo2*((fym2 - fyb) + (fym2 - fya)))
#define BT fym1mfym2 = fym1 - fym2; fym2mfym3 = fym2 - fym3;
#define CROSS1 std::abs(dyo4*(fym1mfym2 + (fym1 - fya))) // see if I can save on computation here
#define CROSS2 std::abs(dyo4*(fym2mfym3 - fym1mfym2))
#define CROSS3 std::abs(dyo4*((fym3 - fyb) - fym2mfym3))
#define CCROSSr1 std::abs(dyo4*(fym1mfym2._real + (fym1._real - fya._real))) // see if I can save on computation here
#define CCROSSr2 std::abs(dyo4*(fym2mfym3._real - fym1mfym2._real))
#define CCROSSr3 std::abs(dyo4*((fym3._real - fyb._real) - fym2mfym3._real))
#define CCROSSi1 std::abs(dyo4*(fym1mfym2._imag + (fym1._imag - fya._imag))) // see if I can save on computation here
#define CCROSSi2 std::abs(dyo4*(fym2mfym3._imag - fym1mfym2._imag))
#define CCROSSi3 std::abs(dyo4*((fym3._imag - fyb._imag) - fym2mfym3._imag))
#define CROSS(num) \
        if constexpr (std::same_as<R, gtd::complex<T>>) { \
            if (CCROSSr1 <= ptol_y && CCROSSr2 <= ptol_y && CCROSSr3 <= ptol_y && \
                CCROSSi1 <= ptol_y && CCROSSi2 <= ptol_y && CCROSSi3 <= ptol_y) {goto divide_##num;} \
        } else { \
            if (CROSS1 <= ptol_y && CROSS2 <= ptol_y && CROSS3 <= ptol_y) {goto divide_##num;} \
        }
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
        auto fy_lam = [&f, &y](const T &x){return f(x, y);};
        using FL = decltype(fy_lam);
        R fyaPfyb;
        R estimate;
        R ires;
        R res = 0;
        bool left = false;
        uint64_t bitmask;
        T absval;
        R fya;
        R fym1;
        R fym2;
        R fym3;
        R fyb;
        R fym1mfym2;
        R fym2mfym3;
        uint64_t *mptr;
#define SIMPDBLQUAD_MOST(par0, par00, par1, par2, par3, lab1, lab2, after) \
        y = ya; \
        fya = simpquad<T, R, FL, false>(fy_lam, gfunc(ya), hfunc(ya), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym1; \
        fym1 = simpquad<T, R, FL, false>(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym2; \
        fym2 = simpquad<T, R, FL, false>(fy_lam, gfunc(ym2), hfunc(ym2), abstol_x, reltol_x, ptol_x, mptr); after \
        y = ym3; \
        fym3 = simpquad<T, R, FL, false>(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr); after \
        y = yb; \
        fyb = simpquad<T, R, FL, false>(fy_lam, gfunc(yb), hfunc(yb), abstol_x, reltol_x, ptol_x, mptr); after \
        while (1) { \
            par1 \
            fyaPfyb = fya + fyb; \
            estimate = dyo6*(fyaPfyb + 4*fym2); \
            ires = dyo12*(fyaPfyb + 4*fym1 + 2*fym2 + 4*fym3); \
            par2 \
            if ((absval = gtd::abs(ires - estimate)) > abstol_y && absval > gtd::abs(reltol_y*ires)) { \
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
                fym1 = simpquad<T, R, FL, false>(fy_lam, gfunc(ym1), hfunc(ym1), abstol_x, reltol_x, ptol_x, mptr);after\
                y = ym3; \
                fym3 = simpquad<T, R, FL, false>(fy_lam, gfunc(ym3), hfunc(ym3), abstol_x, reltol_x, ptol_x, mptr);after\
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
                    fym1 = simpquad<T, R, FL, false>(fy_lam,gfunc(ym1),hfunc(ym1),abstol_x,reltol_x,ptol_x,mptr);after\
                    y = ym3; \
                    fym3 = simpquad<T, R, FL, false>(fy_lam,gfunc(ym3),hfunc(ym3),abstol_x,reltol_x,ptol_x,mptr);after\
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
                        /* par00 */ \
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
                    fym1 = simpquad<T, R, FL, false>(fy_lam,gfunc(ym1),hfunc(ym1),abstol_x,reltol_x,ptol_x,mptr);after\
                    y = ym3; \
                    fym3 = simpquad<T, R, FL, false>(fy_lam,gfunc(ym3),hfunc(ym3),abstol_x,reltol_x,ptol_x,mptr);after\
                    fym3fyb = stack.top(); \
                    fym2 = fym3fyb.fym3; \
                    fyb = fym3fyb.fyb; \
                    left = false; \
                    tree.top() |= bitmask; \
                    par00 \
                } \
            } \
        }
        if (ptol_y >= 0) { // I'm doing this macro insanity to avoid a shit-tonne of code duplication
            const T org_ptol_y = ptol_y;
            if (mdepth_y) {
                // uint64_t _s;
                if (*mdepth_y == NO_MAX_DEPTH) {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                        SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);, EMPTY, BT
                                         CROSS(11),
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_11, other_11, XDEPTH)
                    } else {
                        mptr = nullptr;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8}; // Already provides a depth of 8*64 = 512
                        SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);, EMPTY, BT
                                         CROSS(12),
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
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
                        SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);,
                                         if (stack._size == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_21;}, BT
                                         CROSS(21),
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_21, other_21, XDEPTH)
                    } else {
                        mptr = nullptr;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);,
                                         if (stack._size == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_22;}, BT
                                         CROSS(22),
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_22, other_22, EMPTY)
                    }
                }
            } else {
                if (mdepth_x) {
                    uint64_t maxx_depth = *mdepth_x;
                    const uint64_t org_mxd = *mdepth_x;
                    *mdepth_x = 0;
                    mptr = &maxx_depth;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);, EMPTY, BT
                                     CROSS(31),
                                     EMPTY, divide_31, other_31, XDEPTH)
                } else {
                    mptr = nullptr;
                    gtd::stack<uint64_t> tree{8};
                    SIMPDBLQUAD_MOST(ptol_y /= 4;, ptol_y = org_ptol_y/(1 << 2*stack._size);, EMPTY, BT
                                     CROSS(32),
                                     EMPTY, divide_32, other_32, EMPTY)
                }
            }
        } else {
            if (mdepth_y) {
                // uint64_t _s;
                if (*mdepth_y == NO_MAX_DEPTH) {
                    if (mdepth_x) {
                        uint64_t maxx_depth = *mdepth_x;
                        const uint64_t org_mxd = *mdepth_x;
                        *mdepth_x = 0;
                        mptr = &maxx_depth;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY,
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_41, other_41, XDEPTH)
                    } else {
                        mptr = nullptr;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{8};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY, EMPTY, EMPTY,
                                         if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
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
                                         if (stack._size == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_51;},
                                         EMPTY, if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_51, other_51, XDEPTH)
                    } else {
                        mptr = nullptr;
                        uint64_t max_depth = *mdepth_y;
                        *mdepth_y = 0;
                        gtd::stack<uint64_t> tree{max_depth % 64 ? max_depth/64 + 1 : max_depth/64};
                        SIMPDBLQUAD_MOST(EMPTY, EMPTY,
                                         if (stack._size == max_depth)
                                         {ires = dyo12*(fya + 4*fym1 + 2*fym2 + 4*fym3 + fyb); goto other_52;},
                                         EMPTY, if (stack._size > *mdepth_y) {*mdepth_y = stack._size;},
                                         divide_52, other_52, EMPTY)
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
#undef CCROSSi3
#undef CCROSSi2
#undef CCROSSi1
#undef CCROSSr3
#undef CCROSSr2
#undef CCROSSr1
#undef CROSS1
#undef CROSS2
#undef CROSS3
#undef BT
        if constexpr (prog)
            putchar('\n');
        return res; // for completeness, but would never be reached
    }
}
#endif
