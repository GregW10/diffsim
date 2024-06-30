#include "diffimg.hpp"
#include "../glib/misc/gregparse.hpp"

template <typename T> requires (std::is_floating_point_v<T>)
void gen_bmp(const char *dffr_s, std::string &bmp_path, const char *cmap_str) {
    diff::diffimg<T> img{dffr_s, false};
    diff::colourmap<T> cmap; // grayscale by default
    if (cmap_str) {
        if (auto it = diff::cmaps<T>::all_cmaps.find(cmap_str); it != diff::cmaps<T>::all_cmaps.end())
            cmap = it->second;
        else
            cmap.from_cmap(cmap_str);
    }
    img.gen_bmp(bmp_path, cmap);
}

int convert(gtd::parser &parser) {
    const char *cmap_s = parser.get_arg("--cmap");
    const char *bmp_s = parser.get_arg("-o");
    const char *dffr_s = parser.get_arg(std::regex{R"(^.+$)"});
    if (!dffr_s) {
        fprintf(stderr, "Error: no .dffr file passed.\n");
        return 1;
    }
    std::string path;
    if (bmp_s)
        path = bmp_s;
	try {
        try {
            gen_bmp<long double>(dffr_s, path, cmap_s);
        }
        catch(const diff::invalid_dffr_sizeT&) {
            try {
                gen_bmp<double>(dffr_s, path, cmap_s);
            }
            catch(const diff::invalid_dffr_sizeT&) {
                try {
                    gen_bmp<float>(dffr_s, path, cmap_s);
                }
                catch(const diff::invalid_dffr_sizeT &e) {
                    std::cerr << "Error: an exception occurred. what(): \"" << e.what() << "\"\n";
                    return 1;
                }
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "Error: an exception occurred. what(): \"" << e.what() << "\"\n";
        return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
	gtd::parser parser{argc, argv};
    return convert(parser);
}
