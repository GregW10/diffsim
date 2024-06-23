#include "../utils/diffimg.hpp"

int main() {
    std::cout << diff::cmaps<long double>::grayscale << std::endl;
    diff::cmaps<long double>::grayscale.to_cmap("grayscale_ld.cmap");
    diff::colourmap<long double> bgr{"grayscale_ld.cmap"};
    std::cout << bgr << std::endl;
    diff::colourmap<float> shit;
    return 0;
}
