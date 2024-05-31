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
#ifdef __CUDACC__
    template <numeric T = double>
#else
    template <numeric T = long double>
#endif
    class diff_alloc {
    protected:
        uint64_t w{}; // width of detector (in pixels)
        uint64_t h{}; // height of detector (in pixels)
        uint64_t s{}; // total number of pixels in detector
        uint64_t nb = s*sizeof(T); // total number of bytes allocated
    private:
        gtd::mmapper mapper{s*sizeof(T)}; // gtd::mmapper object to take care of allocation using mmap
    protected:
        T *data = (T*) mapper.get(); // pointer to data on host
#ifdef __CUDACC__
        T *gdat{}; // pointer to data on device
        bool on_gpu = false; // boolean indicating whether stored results are on the GPU
#endif
        bool on_cpu = false; // boolean indicating whether stored results are on the CPU
    public:
        // diff_alloc() = default;
        diff_alloc(uint64_t width, uint64_t height) : w{width}, h{height}, s{width*height} {
#ifdef __CUDACC__
            CUDA_ERROR(cudaMalloc(&gdat, nb));
            CUDA_ERROR(cudaMemset(gdat, 0, nb)); // zeros-out all allocated GPU memory
#endif
            mapper.zero(); // zeros-out all allocated memory
        }
#ifdef __CUDACC__
        bool dev_to_host() {
            if (!on_gpu)
                return false; // no results on GPU, so nothing to copy
            if (on_cpu)
                return true; // already on CPU, so no need to copy
            CUDA_ERROR(cudaMemcpy(data, gdat, nb, cudaMemcpyDeviceToHost));
            on_cpu = true; // now results are on both the CPU and GPU
            return true;
        }
#endif
        uint64_t dtr_width() const noexcept {
            return this->w;
        }
        uint64_t dtr_height() const noexcept {
            return this->h;
        }
        uint64_t dtr_pixels() const noexcept {
            return this->s;
        }
        uint64_t nbytes() const noexcept {
            return this->nb;
        }
        uint64_t to_dttr(const char *path) {
            if (!path || !*path)
                return 0;
#ifdef __CUDACC__
            if (!dev_to_host()) { // takes care of copying results to CPU if on GPU, and if not, returns false
                return 0;
            }
#else
            if (!on_cpu)
                return 0; // because results have not yet been calculated
#endif
            int fd;
            if ((fd = open(path, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR)) == -1)
                return 0;
            uint64_t bwritten = 0;
            char hdr[4] = {'D', 'T', 'T', 'R'};
            uint64_t sizeT = sizeof(T);
            if (!gtd::write_all(fd, hdr, 4))
                goto end;
            bwritten += 4;
            if (!gtd::write_all(fd, &this->w, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, &this->h, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, &sizeT, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, this->data, this->nb))
                goto end;
            bwritten += this->nb;
            end:
            close(fd);
            return bwritten;
        }
        uint64_t from_dttr(const char *path) {
            if (!path || !*path)
                return 0;
            struct stat buff{};
            if (stat(path, &buff) == -1)
                return 0;
            int fd;
            if ((fd = open(path, O_RDONLY)) == -1)
                return 0;
            uint64_t bwritten = 0;
            char hdr[4] = {'D', 'T', 'T', 'R'};
            uint64_t sizeT = sizeof(T);
            if (!gtd::write_all(fd, hdr, 4))
                goto end;
            bwritten += 4;
            if (!gtd::write_all(fd, &this->w, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, &this->h, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, &sizeT, sizeof(uint64_t)))
                goto end;
            bwritten += sizeof(uint64_t);
            if (!gtd::write_all(fd, this->data, this->nb))
                goto end;
            bwritten += this->nb;
            end:
            close(fd);
            return bwritten;
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
