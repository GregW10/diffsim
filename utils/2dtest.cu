#include "diffbase.hpp"

#define DBL double

#define DFORMAT "lf"

#define NUM_BLOCKS 2048

class Functor2D {
    mutable uint64_t times_called = 0;
public:
    HOST_DEVICE DBL operator()(DBL x, DBL y) const noexcept {
        ++times_called;
        return std::sin(x)*std::cos(y);
    }
    __host__ __device__ uint64_t reset() const noexcept {
        uint64_t val = times_called;
        times_called = 0;
        return val;
    }
    __host__ __device__ uint64_t ncalls() const noexcept {
        return times_called;
    }
};

__host__ __device__ DBL gfunc(DBL y) {
    y -= 1;
    //printf("gfunc: %lf\n", 1 - std::sqrt(1 - y*y));
    return 1 - std::sqrt(1 - y*y);
}

__host__ __device__ DBL hfunc(DBL y) {
    y -= 1;
    //printf("hfunc: %lf\n", 1 + std::sqrt(1 - y*y));
    return 1 + std::sqrt(1 - y*y);
}

__global__ void kernel_main(DBL *res, uint64_t *ncalls, uint64_t *mdepth_x, uint64_t *mdepth_y) {
    printf("inside kernel\n");
    Functor2D f;
    DBL ya = 0;
    DBL yb = 2;
    DBL abstol_y = 0.0625l/1024000.0;
    DBL reltol_y = 0.0625l/1024000.0;
    DBL ptol_y = 0.00000000000000001;
    DBL abstol_x = 0.0625l/102400000.0;
    DBL reltol_x = 0.0625l/102400000.0;
    DBL ptol_x = 0.00000000000000001;
    *(res + blockIdx.x) = diff::simpdblquad(f, ya, yb, gfunc, hfunc, abstol_y, reltol_y, ptol_y,
                                            mdepth_y + blockIdx.x,
                                abstol_x, reltol_x, ptol_x, mdepth_x + blockIdx.x);
    *(ncalls + blockIdx.x) = f.ncalls();
}

void print_error(cudaError_t err) {
    if (err != cudaSuccess) {
        fprintf(stderr, "Error: %s\n", cudaGetErrorString(err));
        exit(1);
    }
}

int main(int argc, char **argv) {
    uint64_t mdepth_x = 36;
    uint64_t mdepth_y = 36;
    if (argc == 2)
        mdepth_x = mdepth_y = diff::misc::to_uint(*(argv + 1));
    else if (argc == 3) {
        mdepth_x = diff::misc::to_uint(*(argv + 1));
        mdepth_y = diff::misc::to_uint(*(argv + 2));
    }
    else if (argc > 3)
        return 1;
    uint64_t mdx_host[NUM_BLOCKS];
    uint64_t mdy_host[NUM_BLOCKS];
    uint64_t *xptr = mdx_host;
    uint64_t *yptr = mdy_host;
    unsigned int i = NUM_BLOCKS;
    while (i --> 0)
        *xptr++ = mdepth_x;
    i = NUM_BLOCKS;
    while (i --> 0)
        *yptr++ = mdepth_y;
    DBL *res_device;
    print_error(cudaMalloc(&res_device, NUM_BLOCKS*sizeof(DBL)));
    uint64_t *ncalls_device;
    cudaMalloc(&ncalls_device, NUM_BLOCKS*sizeof(uint64_t));
    uint64_t *mdx_device;
    cudaMalloc(&mdx_device, NUM_BLOCKS*sizeof(uint64_t));
    uint64_t *mdy_device;
    cudaMalloc(&mdy_device, NUM_BLOCKS*sizeof(uint64_t));
    cudaMemcpy(mdx_device, mdx_host, NUM_BLOCKS*sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(mdy_device, mdy_host, NUM_BLOCKS*sizeof(uint64_t), cudaMemcpyHostToDevice);
    kernel_main<<<NUM_BLOCKS,1>>>(res_device, ncalls_device, mdx_device, mdy_device);
    cudaDeviceSynchronize();
    DBL res_host[NUM_BLOCKS];
    uint64_t ncalls_host[NUM_BLOCKS];
    cudaMemcpy(res_host, res_device, NUM_BLOCKS*sizeof(DBL), cudaMemcpyDeviceToHost);
    cudaMemcpy(ncalls_host, ncalls_device, NUM_BLOCKS*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(mdx_host, mdx_device, NUM_BLOCKS*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(mdy_host, mdy_device, NUM_BLOCKS*sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaFree(res_device);
    cudaFree(ncalls_device);
    cudaFree(mdx_device);
    cudaFree(mdy_device);
    i = NUM_BLOCKS;
    DBL *dptr = res_host;
    uint64_t *nptr = ncalls_host;
    xptr = mdx_host;
    yptr = mdy_host;
    while (i --> 0) {
        printf("Result: %.20" DFORMAT "\n", *dptr++);
        printf("Num func. calls: %" PRIu64 "\n", *nptr++);
        printf("Max. depth in y: %" PRIu64 "\n", *yptr);
        printf("Max. depth in x: %" PRIu64 "\n", *xptr);
    }
    return 0;
}
