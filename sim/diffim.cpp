#include "../utils/diffimg.hpp"
#include "../glib/misc/gregparse.hpp"

template <typename T> requires (std::is_floating_point_v<T>)
void start_sim(gtd::parser &parser) {

}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};

    return 0;
}
