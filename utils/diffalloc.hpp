#ifndef DIFFALLOC_HPP
#define DIFFALLOC_HPP

#include "diffbase.hpp"
#include "../glib/misc/gregmisc.hpp"
#include "../glib/misc/gregmmapper.hpp"
#include <fcntl.h>

namespace diff {
#ifdef __CUDACC__
    class cuda_error : public std::runtime_error {
        char *msg{};
    public:
        explicit cuda_error(int line, const char *file, const char *error_s) {
            msg = new char[40 + gtd::strlen_c(file) + gtd::strlen_c(error_s)];
            gtd::strcpy_c(msg, "Error \"");
            gtd::strcat_c(msg, error_s);
            gtd::strcat_c(msg, "\" occurred in file \"");
            gtd::strcat_c(msg, file);
            gtd::strcat_c(msg, "\" on line ");
            gtd::to_string(line, msg);
            gtd::strcat_c(msg, ".\n");
        }
        const char *what() const noexcept override {
            return this->msg;
        }
        ~cuda_error() {
            delete [] msg;
        }
    };
#define CUDA_ERROR(func_call) \
    cudaError_t err; \
    if ((err = func_call) != cudaSuccess) { \
        throw cuda_error{__LINE__, __FILE__, cudaGetErrorString(err)}; \
    }
#endif
    template <numeric T>
    class diff_alloc {
    protected:
        uint64_t w{}; // width of detector (in pixels)
        uint64_t h{}; // height of detector (in pixels)
        uint64_t s{}; // total number of pixels in detector
    private:
        gtd::mmapper mapper{s*sizeof(T)}; // gtd::mmapper object to take care of allocation using mmap
    protected:
        T *data = mapper.get(); // pointer to data on host
#ifdef __CUDACC__
        T *gdat{}; // pointer to data on device
        bool on_gpu = false; // boolean indicating whether stored results are on the GPU
#endif
        bool on_cpu = false; // boolean indicating whether stored results are on the CPU
    public:
        diff_alloc() = default;
        diff_alloc(uint64_t width, uint64_t height) : w{width}, h{height}, s{width*height} {
#ifdef __CUDACC__
            CUDA_ERROR(cudaMalloc(&gdat, sizeof(T)*s));
            CUDA_ERROR(cudaMemset(gdat, 0, sizeof(T)*s)); // zeros-out all allocated GPU memory
#endif
            mapper.zero(); // zeros-out all allocated memory
        }

#ifdef __CUDACC__
        bool dev_to_host() {
            if (!on_gpu)
                return false;
            CUDA_ERROR(cudaMemcpy(data, gdat, sizeof(T)*s, cudaMemcpyDeviceToHost));
            on_gpu = false; // since stored results are no longer ONLY on the GPU
            return true;
        }
#endif
        bool to_dttr(const char *path) {
            if (!path || !*path)
                return false;
#ifdef __CUDACC__
            if (on_gpu) {

            }
#endif
            int fd;
            if ((fd = open(path, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR)) == -1)
                return false;
            char hdr[4] = {'D', 'T', 'T', 'R'};
            if (!gtd::write_all(fd, hdr, 4))
                return false;

        }
        ~diff_alloc() { // gtd::mmapper object takes care of the deallocation in its destructor
            // delete [] data;
#ifdef __CUDACC__
            cudaFree(gdat); // no error handling, as I should not be throwing inside a destructor!
#endif
        }
    };
}
#endif
