#ifndef DIFFBASE_HPP
#define DIFFBASE_HPP

#ifdef __CUDACC__
#warning "Compiling with nvcc.\n"
#define HOST_DEVICE __host__ __device__
//#define HOST __host__
#define DEVICE __device__
#define GLOBAL __global__
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
    HOST_DEVICE inline T abs(const T &num) {
        return num >= 0 ? num : -num;
    }
    template <numeric T>
    HOST_DEVICE T integrate_trap(T (*f)(const T&), const T &a, const T &b, uint64_t num = 1'000'000) {
        if (!f)
            return 0;
        T dx = (b - a)/num;
        T dx_o2 = dx/2;
        uint64_t i = 0;
        uint64_t j = 1;
        T res = 0;
        while (i < num)
            res += dx_o2*(f(a + i++*dx) + f(a + j++*dx));
        return res;
    }
    template <numeric T>
    HOST_DEVICE T trapquad_recurse(T (*f)(const T&), const T &a, const T &b, const T &tol, const T &f_a, const T &f_b) {
        T estimate = ((b - a)/2)*(f_a + f_b);
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_a + 2*f_m + f_b);
        if (abs(res_or_tol - estimate) > tol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T>
    HOST_DEVICE T integrate_trapquad(T (*f)(const T&), const T &a, const T &b, const T &tol = 0.0625l/1024.0l) {
        if (!f)
            return 0;
        T f_a = f(a);
        T f_b = f(b);
        T estimate = ((b - a)/2)*(f_a + f_b);
        T midpoint = (a + b)/2;
        T f_m = f(midpoint);
        T res_or_tol = ((b - midpoint)/2)*(f_a + 2*f_m + f_b);
        if (abs(res_or_tol - estimate) > tol) {
            res_or_tol = tol/2;
            return trapquad_recurse(f, a, midpoint, res_or_tol, f_a, f_m) +
                   trapquad_recurse(f, midpoint, b, res_or_tol, f_m, f_b);
        }
        return res_or_tol;
    }
    template <numeric T>
    HOST_DEVICE T integrate_simpson(T (*f)(const T&), const T &a, const T &b, uint64_t num = 1'000'000) {
        if (!f)
            return 0;
        T dx = (b - a)/num;
        // num *= 2;
        //T dx_o2 = dx/2;//(b - a)/(num);
        uint64_t i = 1;
        T offset;
        T a_p_dxo2 = a + dx/2;//dx_o2;
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
