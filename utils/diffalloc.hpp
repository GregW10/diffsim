#ifndef DIFFALLOC_HPP
#define DIFFALLOC_HPP

#include "diffbase.hpp"

namespace diff {
#ifdef __CUDACC__
    #define CUDA_ERROR(func_call) \
    cudaError_t err; \
    if ((err = func_call) != cudaSuccess) { \
        fprintf(stderr, "CUDA error on line %d in file \"%s\": %s\n", __LINE__, __FILE__, cudaGetErrorString(err)); \
        exit(1); \
    }
#endif
    template <numeric T>
    class diff_alloc {
        uint64_t w{}; // width of detector (in pixels)
        uint64_t h{}; // height of detector (in pixels)
        uint64_t s{}; // total number of pixels in detector
        T *data{}; // pointer to data on host
#ifdef __CUDACC__
        T *gdat{}; // pointer to data on device
#endif
    public:
        diff_alloc() = default;
        diff_alloc(uint64_t width, uint64_t height) : w{width}, h{height}, s{width*height}, data{new T[s]} {
#ifdef __CUDACC__
            CUDA_ERROR(cudaMalloc(&gdat, sizeof(T)*s));
#endif
        }
        ~diff_alloc() {
            delete [] data;
#ifdef __CUDACC__

#endif
        }
    };
}
#endif
