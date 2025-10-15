#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "secret.hpp"
#include "image.hpp"

namespace fs = std::filesystem;

int main() {
    const std::string input_path = "images/input.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
    size_t secret_size = secret_pixels.size();
    std::cout << "Loaded secret image '" << input_path << "' (" << width << "x" << height << "), bytes: " << secret_size << "\n";

    const unsigned n = 3;
    const unsigned k = 3;
    const unsigned kn = 1;

    const std::string out_base = "results/shamir256jpeg_x";
    fs::create_directories(out_base);
    fs::create_directories(out_base + "/shadows");
    fs::create_directories(out_base + "/shadows_jpeg");
    fs::create_directories(out_base + "/reconstructions");
    fs::create_directories(out_base + "/diff_maps");

    auto field = ss::BinaryField<std::uint8_t>(8);

    std::vector<std::vector<std::uint8_t>> x_sets = {
        {1, 2, 3},
        {1, 255, 254},
        {1, 128, 129},
        {1, 2, 128},
        {1, 3, 5},
        {85, 170, 255}
    };

    const int jpeg_quality = 100;

    for (size_t set_idx = 0; set_idx < x_sets.size(); ++set_idx) {
        const auto& x_values = x_sets[set_idx];
        std::cout << "\n=== Testing x_k set " << (set_idx + 1) << ": {";
        for (size_t i = 0; i < x_values.size(); ++i) {
            std::cout << static_cast<int>(x_values[i]);
            if (i < x_values.size() - 1) std::cout << ", ";
        }
        std::cout << "} ===\n";

        std::vector<std::uint8_t> field_vec = secret_pixels;
        auto shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn, x_values);

        if (shares.size() != n) {
            std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
            continue;
        }

        for (size_t i = 0; i < shares.size(); ++i) {
            std::string outp = out_base + "/shadows/set" + std::to_string(set_idx+1) +
                              "_shadow_x" + std::to_string(static_cast<int>(x_values[i])) + ".png";
            ss::saveGrayscalePNG(outp, shares[i], width, height);
        }

        for (unsigned share_to_noise = 0; share_to_noise < n; ++share_to_noise) {
            std::cout << "\n--- Noising share " << (share_to_noise + 1) << " (x="
                      << static_cast<int>(x_values[share_to_noise]) << ") ---\n";

            std::vector<std::vector<std::uint8_t>> comp_shares = shares;

            std::string jpeg_path = out_base + "/shadows_jpeg/set" + std::to_string(set_idx+1) +
                                   "_noised_share" + std::to_string(share_to_noise+1) +
                                   "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".jpg";

            if (!stbi_write_jpg(jpeg_path.c_str(), width, height, 1, shares[share_to_noise].data(), jpeg_quality)) {
                std::cerr << "Warning: failed to write JPEG for share " << (share_to_noise+1) << " to '" << jpeg_path << "'\n";
                continue;
            }

            int w2=0, h2=0, channels=0;
            unsigned char *data = stbi_load(jpeg_path.c_str(), &w2, &h2, &channels, 1);
            if (!data) {
                std::cerr << "Warning: failed to reload JPEG '" << jpeg_path << "': " << stbi_failure_reason() << "\n";
                continue;
            }

            if (w2 != width || h2 != height) {
                std::cerr << "Warning: reloaded JPEG size mismatch for '" << jpeg_path << "' (" << w2 << "x" << h2 << ") expected (" << width << "x" << height << ")\n";
            }
            comp_shares[share_to_noise].assign(data, data + (static_cast<size_t>(w2) * static_cast<size_t>(h2)));
            stbi_image_free(data);

            std::cout << "Share " << (share_to_noise+1) << " (x=" << static_cast<int>(x_values[share_to_noise])
                      << ") orig vs JPEG | PSNR = " << ss::computePSNR(shares[share_to_noise], comp_shares[share_to_noise])
                      << " dB | NPCR = " << ss::computeNPCR(shares[share_to_noise], comp_shares[share_to_noise])
                      << "% | UACI = " << ss::computeUACI(shares[share_to_noise], comp_shares[share_to_noise]) << "%\n";

            std::string diff_path = out_base + "/diff_maps/diff_set" + std::to_string(set_idx+1) +
                                   "_share" + std::to_string(share_to_noise+1) +
                                   "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";

            std::vector<std::uint8_t> diff_map = ss::generateDiffMap(shares[share_to_noise], comp_shares[share_to_noise], width, height);
            ss::saveGrayscalePNG(diff_path, diff_map, width, height);

            // Choose k shares, prefer non-noised ones
            std::vector<unsigned> indices;
            indices.reserve(k);
            for (unsigned i = 0; i < comp_shares.size() && indices.size() < k; ++i) {
                if (i == share_to_noise) continue;
                indices.push_back(i);
            }
            if (indices.size() < k) {
                for (unsigned i = 0; i < comp_shares.size() && indices.size() < k; ++i) {
                    if (std::find(indices.begin(), indices.end(), i) == indices.end())
                        indices.push_back(i);
                }
            }

            auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints<std::uint8_t>(indices, comp_shares, x_values);
            auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));

            std::string rec_path = out_base + "/reconstructions/reconstruction_set" + std::to_string(set_idx+1) +
                                  "_noised_share" + std::to_string(share_to_noise+1) +
                                  "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";
            ss::saveGrayscalePNG(rec_path, rec, width, height);

            double psnr = ss::computePSNR(field_vec, rec);
            double npcr = ss::computeNPCR(field_vec, rec);
            double uaci = ss::computeUACI(field_vec, rec);

            std::cout << "Reconstruction with x_set {"
                      << static_cast<int>(x_values[0]) << ", "
                      << static_cast<int>(x_values[1]) << ", "
                      << static_cast<int>(x_values[2])
                      << "}, noised share " << (share_to_noise+1)
                      << " | PSNR: " << psnr << " dB | NPCR: " << npcr
                      << "% | UACI: " << uaci << "%\n";

            std::string rec_diff = out_base + "/diff_maps/diff_reconstruction_set" + std::to_string(set_idx+1) +
                                  "_noised_share" + std::to_string(share_to_noise+1) +
                                  "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";

            std::vector<std::uint8_t> rec_map = ss::generateDiffMap(field_vec, rec, width, height);
            ss::saveGrayscalePNG(rec_diff, rec_map, width, height);
        }
    }

    ss::saveHistogram(out_base + "/histogram.dat", ss::computeHistogram(secret_pixels));

    std::cout << "\n=== All tests completed. Total: " << x_sets.size() << " sets × 3 shares = "
              << (x_sets.size() * 3) << " test cases ===\n";

    return 0;
}