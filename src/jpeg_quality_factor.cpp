#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <map>

#include "ss/image.hpp"

namespace fs = std::filesystem;

int main() {
    std::vector<std::string> input_images = {
        "images/baboon.png",
        "images/frog.png",
        "images/louvre.png",
        "images/peppers.png",
        "images/stop.png"
    };

    const int min_quality = 1;
    const int max_quality = 100;

    const std::string npcr_file = "results/jpeg_quality/npcr.dat";
    fs::create_directories(fs::path(npcr_file).parent_path());

    const std::string hist_dir = "results/jpeg_quality/diff_hist/";
    fs::create_directories(hist_dir);

    std::ofstream npcr_ofs(npcr_file);
    if (!npcr_ofs.is_open()) {
        std::cerr << "Failed to open " << npcr_file << " for writing\n";
        return 1;
    }

    struct ImageData {
        std::vector<std::uint8_t> pixels;
        int width;
        int height;
    };

    std::vector<ImageData> images;
    for (const auto &img_path : input_images) {
        auto [pixels, width, height] = ss::readGrayscalePNG(img_path);
        images.push_back(ImageData{pixels, width, height});
    }

    for (int quality = max_quality; quality >= min_quality; --quality) {
        double sum_npcr = 0.0;
        std::map<int, int> diff_hist;
        size_t total_pixels = 0;

        for (const auto &img : images) {
            auto jpeg_pixels = ss::jpegify(img.pixels, img.width, img.height, quality);

            double npcr = ss::computeNPCR(img.pixels, jpeg_pixels);
            sum_npcr += npcr;

            if (quality % 5 == 0) {
                for (size_t i = 0; i < img.pixels.size(); ++i) {
                    int diff = static_cast<int>(jpeg_pixels[i]) - static_cast<int>(img.pixels[i]);
                    diff_hist[diff]++;
                }
                total_pixels += img.pixels.size();
            }
        }

        double avg_npcr = sum_npcr / images.size();
        npcr_ofs << quality << " " << avg_npcr << "\n";
        std::cout << quality << " " << avg_npcr << "\n";

        if (quality % 5 == 0) {
            std::string hist_file = hist_dir + "diff_hist" + std::to_string(quality) + ".dat";
            std::ofstream hist_ofs(hist_file);
            if (!hist_ofs.is_open()) {
                continue;
            }

            for (const auto &[diff, count] : diff_hist) {
                double percentage = 100.0 * count / total_pixels;
                hist_ofs << diff << " " << percentage << "\n";
            }

            hist_ofs.close();
        }
    }

    npcr_ofs.close();
    return 0;
}