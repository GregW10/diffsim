#include "../glib/misc/gregstack.hpp"
#include <iostream>

__global__ void kernel(size_t size) {
    gtd::stack<int> stack;
    printf("Stack created...\n");
    int i = 0;
    while (i < 129)
        stack.push(i++);
    printf("Capacity: %" PRIu64 "\n", stack.capacity());
    while (stack)
        printf("Int: %d\n", stack.pop_r());
}

int main() {
    dim3 blocks{16, 16};
    dim3 grid{1, 1};
    size_t size = 1'000;
    kernel<<<grid,blocks>>>(size);
    cudaDeviceSynchronize();
    return 0;
}
