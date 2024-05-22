#include "../glib/misc/gregstack.hpp"
#include <iostream>

int main() {
    gtd::stack<int> stack;
    std::cout << "Stack created.\n" << std::endl;
    int i = 0;
    while (i < 129)
        stack.push(i++);
    std::cout << "Capacity: " << stack.capacity() << std::endl;
    while (stack)
        std::cout << stack.pop_r() << std::endl;
    return 0;
}
