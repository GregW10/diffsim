#include "../glib/misc/gregcomplex.hpp"

int main() {
    gtd::complex<float> c{1, 1};
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.real(-1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.imag(-1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c.real(1);
    std::cout << c << std::endl;
    std::cout << +gtd::rad_to_deg(c.arg()) << std::endl;
    c = c*5;
    std::cout << c << std::endl;
    c = 5*c;
    std::cout << c << std::endl;
    return 0;
}
