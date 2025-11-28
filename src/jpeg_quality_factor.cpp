#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <numeric>
#include <fstream>
#include <map>
#include <algorithm>
#include <random>
#include <iomanip>

#include "ss/image.hpp"
#include "ss/secret.hpp"

namespace fs = std::filesystem;

int main() {
    using Storage = std::uint8_t;

    const int min_quality = 1;
    const int max_quality = 100;
    const unsigned shk = 10;
    const unsigned kn = 1;
    const unsigned max_images = 20;

    const std::string input_dir = "images/minibows2/";
    const std::string out_dir = "results/jpeg_quality/";
    const std::string hist_dir = out_dir + "diff_hist/";

    fs::create_directories(out_dir);
    fs::create_directories(hist_dir);

    std::vector<std::string> input_images;
    unsigned img_count = 0;
    for (auto const &entry : fs::directory_iterator(input_dir)) {
        if (img_count >= max_images) break;
        if (!entry.is_regular_file()) continue;
        input_images.push_back(entry.path().string());
        ++img_count;
    }

    if (input_images.empty()) {
        std::cerr << "No PNG files found in " << input_dir << "\n";
        return 1;
    }

    ss::BinaryField<Storage> F = ss::BinaryField<Storage>(8);

    struct ImageData {
        std::vector<std::uint8_t> pixels;
        int width;
        int height;
        std::string path;
    };

    std::vector<ImageData> images;
    images.reserve(input_images.size());
    for (const auto &img_path : input_images) {
        auto [pixels, width, height] = ss::readGrayscale(img_path);
        images.push_back(ImageData{pixels, width, height, img_path});
    }

    std::ofstream npcr_ofs(out_dir + "npcr.dat");
    npcr_ofs << std::fixed << std::setprecision(6);

    for (int quality = max_quality; quality >= min_quality; --quality) {
        double sum_npcr = 0.0;
        std::uint64_t count_npcr_items = 0;
        std::map<int, std::uint64_t> diff_hist;
        std::uint64_t total_pixels = 0;

        for (const auto &img : images) {
            std::vector<Storage> data;
            data.reserve(img.pixels.size());
            for (auto v : img.pixels) data.push_back(static_cast<Storage>(v));

            unsigned k = shk;
            unsigned n = k;

            std::vector<std::vector<Storage>> shares;
            shares = ss::getShares(data, k, n, F, kn);
            
            // For every share: jpegify and compute NPCR
            for (const auto &share : shares) {
                std::vector<std::uint8_t> share_bytes;
                share_bytes.reserve(share.size());
                for (auto s : share) share_bytes.push_back(static_cast<std::uint8_t>(s));

                auto jpeg_share = ss::jpegify(share_bytes, img.width, img.height, quality);

                double npcr = ss::computeNPCR(share_bytes, jpeg_share);
                sum_npcr += npcr;
                ++count_npcr_items;

                if (quality % 5 == 0) {
                    for (size_t i = 0; i < share_bytes.size(); ++i) {
                        int diff = static_cast<int>(jpeg_share[i]) - static_cast<int>(share_bytes[i]);
                        diff_hist[diff] += 1;
                    }
                    total_pixels += static_cast<std::uint64_t>(share_bytes.size());
                }
            }
        }

        double avg_npcr = (count_npcr_items > 0) ? (sum_npcr / static_cast<double>(count_npcr_items)) : 0.0;
        npcr_ofs << quality << " " << avg_npcr << "\n";
        std::cout << "quality=" << quality << " avg_npcr=" << avg_npcr << "\n";

        if (quality % 5 == 0) {
            std::string hist_file = hist_dir + "diff_hist" + std::to_string(quality) + ".dat";
            std::ofstream hist_ofs(hist_file);
            if (!hist_ofs.is_open()) {
                std::cerr << "Failed to open " << hist_file << " for writing\n";
            } else {
                if (diff_hist.empty()) {
                    hist_ofs << "-1 0\n0 0\n1 0\n";
                } else {
                    int min_diff = diff_hist.begin()->first;
                    int max_diff = diff_hist.rbegin()->first;
                    int left_pad = min_diff - 1;
                    int right_pad = max_diff + 1;

                    for (int d = left_pad; d <= right_pad; ++d) {
                        std::uint64_t cnt = 0;
                        auto it = diff_hist.find(d);
                        if (it != diff_hist.end()) cnt = it->second;
                        double percentage = 0.0;
                        if (total_pixels > 0) percentage = 100.0 * static_cast<double>(cnt) / static_cast<double>(total_pixels);
                        hist_ofs << d << " " << percentage << "\n";
                    }
                }
                hist_ofs.close();
            }
        }
    }

    npcr_ofs.close();

    return 0;
}