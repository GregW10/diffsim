#ifndef DIFFALLOC_HPP
#define DIFFALLOC_HPP

#include "diffbase.hpp"
#include "../glib/misc/gregmisc.hpp"
#include "../glib/misc/gregmmapper.hpp"
#include <fcntl.h>

namespace diff {
#ifdef __CUDACC__
    class cuda_error : public std::exception {
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
    { cudaError_t err; \
    if ((err = func_call) != cudaSuccess) { \
        fprintf(stderr, "Error: %s\n", cudaGetErrorString(err)); \
        throw diff::cuda_error{__LINE__, __FILE__, cudaGetErrorString(err)}; \
    } }
#endif
    struct coord {
        uint64_t x;
        uint64_t y;
    };
#ifdef __CUDACC__
    template <gtd::numeric T = double>
#else
    template <gtd::numeric T = long double>
#endif
    class diffalloc {
    protected:
        uint64_t nw{}; // width of detector (in pixels)
        uint64_t nh{}; // height of detector (in pixels)
        uint64_t np{}; // total number of pixels in detector
        uint64_t nb = np*sizeof(T); // total number of bytes allocated
        gtd::mmapper mapper{np*sizeof(T)}; // gtd::mmapper object to take care of allocation using mmap
        T *data = (T*) mapper.get(); // pointer to data on host
#ifdef __CUDACC__
        T *gdat{}; // pointer to data on device
        // bool on_gpu = false; // boolean indicating whether stored results are on the GPU
#endif
        // bool on_cpu = false; // boolean indicating whether stored results are on the CPU
        HOST_DEVICE uint64_t pix_offset(uint64_t x, uint64_t y) { // get offset into array from pixel x- and y-coords
            return y*nw + x;
        }
        HOST_DEVICE coord pix_coords(uint64_t offset) { // get x- and y-coordinates of pixel at given offset
            return {offset % nw, offset/nw};
        }
    public:
        diffalloc() = default;
        diffalloc(uint64_t width, uint64_t height) : nw{width}, nh{height}, np{width*height} {
#ifdef __CUDACC__
            CUDA_ERROR(cudaMalloc(&gdat, nb));
            CUDA_ERROR(cudaMemset(gdat, 0, nb)); // zeros-out all allocated GPU memory
#endif
            mapper.zero(); // zeros-out all allocated memory
        }
        explicit diffalloc(const char *dttr_path) {
            this->from_dttr(dttr_path);
        }
#ifdef __CUDACC__
        /* bool dev_to_host() {
            if (!on_gpu)
                return false; // no results on GPU, so nothing to copy
            if (on_cpu)
                return true; // already on CPU, so no need to copy
            CUDA_ERROR(cudaMemcpy(data, gdat, nb, cudaMemcpyDeviceToHost));
            on_cpu = true; // now results are on both the CPU and GPU
            return true;
        } */
#endif
        uint64_t dpwidth() const noexcept {
            return this->nw;
        }
        uint64_t dpheight() const noexcept {
            return this->nh;
        }
        uint64_t dpixels() const noexcept {
            return this->np;
        }
        uint64_t dnbytes() const noexcept {
            return this->nb;
        }
        uint64_t dttr_fsize() const noexcept {
            return 28 + this->nb;
        }
        uint64_t zmem(bool just_data = true) noexcept {
            if (!just_data)
                return this->mapper.zero();
            uint64_t counter = this->nb;
            T *ptr = this->data;
            while (counter --> 0)
                *ptr++ = 0;
            return this->nb;
        }
        uint64_t to_dttr(const char *path) {
            if (!path || !*path)
                return 0;
#ifdef __CUDACC__
            // if (!dev_to_host()) // takes care of copying results to CPU if on GPU, and if not, returns false
                // return 0;
            // dev_to_host(); // copy results to CPU if on GPU
#else
            // if (on_cpu)
                // return 0; // because results have not yet been calculated
#endif
            int fd;
            if ((fd = open(path, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR)) == -1)
                return 0;
            uint64_t tot_written = 0;
            uint64_t bwritten;
            char hdr[4] = {'D', 'T', 'T', 'R'};
            uint64_t sizeT = sizeof(T);
            if ((bwritten = gtd::write_all(fd, hdr, 4)) != 4) {
                tot_written += bwritten;
                goto end;
            }
            tot_written += 4;
            if ((bwritten = gtd::write_all(fd, &this->nw, sizeof(uint64_t))) != sizeof(uint64_t)) {
                tot_written += bwritten;
                goto end;
            }
            tot_written += sizeof(uint64_t);
            if ((bwritten = gtd::write_all(fd, &this->nh, sizeof(uint64_t))) != sizeof(uint64_t)) {
                tot_written += bwritten;
                goto end;
            }
            tot_written += sizeof(uint64_t);
            if ((bwritten = gtd::write_all(fd, &sizeT, sizeof(uint64_t))) != sizeof(uint64_t)) {
                tot_written += bwritten;
                goto end;
            }
            tot_written += sizeof(uint64_t);
            if ((bwritten = gtd::write_all(fd, this->data, this->nb)) != this->nb) {
                tot_written += bwritten;
                goto end;
            }
            tot_written += this->nb;
            end:
            close(fd);
            return tot_written;
        }
        uint64_t from_dttr(const char *path) {
            if (!path || !*path)
                return 0;
            struct stat buff{};
            if (stat(path, &buff) == -1)
                return 0;
            if (buff.st_size < 28)
                return 0;
            int fd;
            if ((fd = open(path, O_RDONLY)) == -1)
                return 0;
            uint64_t tot_read = 0;
            uint64_t bread;
            char info[28];
            if ((bread = gtd::read_all(fd, info, 28)) != 28) {
                close(fd);
                return bread;
            }
            tot_read += bread;
            if (info[0] != 'D' || info[1] != 'T' || info[2] != 'T' || info[3] != 'R') {
                close(fd);
                return tot_read;
            }
            if (*((uint64_t*) (info + 20)) != sizeof(T)) {
                close(fd);
                return tot_read;
            }
            this->nw = *((uint64_t*) (info + 4));
            this->nh = *((uint64_t*) (info + 12));
            uint64_t old_np = this->np;
            this->np = this->nw*this->nh;
            this->nb = this->np*sizeof(T);
            if (buff.st_size != 28 + this->nb) {
                close(fd);
                return tot_read;
            }
            if (this->np != old_np) {
                // delete [] data;
                // data = new T[this->s];
                mapper.reset(this->nb);
                data = (T*) mapper.get();
            }
            if ((bread = gtd::read_all(fd, data, this->nb)) != this->nb) {
                close(fd);
                return tot_read + bread;
            }
            close(fd);
            // on_cpu = true;
#ifdef __CUDACC__
            // on_gpu = false;
#endif
            return tot_read;
        }
        ~diffalloc() { // gtd::mmapper object takes care of the deallocation in its destructor
            // delete [] data;
#ifdef __CUDACC__
            cudaFree(gdat); // no error handling, as I should not be throwing inside a destructor!
#endif
        }
    };
}
#endif
