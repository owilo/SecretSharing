#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <limits>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "secret.hpp"
#include "image.hpp"

namespace fs = std::filesystem;

bool saveGrayPNG(const std::string &path, const std::vector<std::uint8_t> &pixels, int w, int h) {
    if ((static_cast<std::size_t>(w) * h) != pixels.size()) {
        std::cerr << "Warning: pixel buffer size (" << pixels.size() << ") != w*h (" << (w*h) << ").\n";
        std::vector<std::uint8_t> tmp((size_t)w * (size_t)h, 0);
        size_t copy = std::min(tmp.size(), pixels.size());
        std::copy(pixels.begin(), pixels.begin() + copy, tmp.begin());
        return stbi_write_png(path.c_str(), w, h, 1, tmp.data(), w) != 0;
    }
    return stbi_write_png(path.c_str(), w, h, 1, pixels.data(), w) != 0;
}

bool readGrayPNG(const std::string &path, std::vector<std::uint8_t> &out_pixels, int &out_w, int &out_h) {
    int w,h,channels;
    unsigned char *data = stbi_load(path.c_str(), &w, &h, &channels, 1);
    if (!data) {
        std::cerr << "Failed to load image '" << path << "': " << stbi_failure_reason() << "\n";
        return false;
    }
    out_w = w; out_h = h;
    out_pixels.assign(data, data + (w * h));
    stbi_image_free(data);
    return true;
}

void selectSharesAndXs(const std::vector<std::vector<std::uint8_t>> &all_shares,
                       const std::vector<unsigned> &indices0based,
                       std::vector<std::vector<std::uint8_t>> &sel_shares,
                       std::vector<std::uint8_t> &xs)
{
    sel_shares.clear();
    xs.clear();
    for (unsigned idx : indices0based) {
        if (idx >= all_shares.size()) {
            std::cerr << "Index out of range when selecting shares: " << idx << "\n";
            continue;
        }
        sel_shares.push_back(all_shares[idx]);
        xs.push_back(static_cast<std::uint8_t>(idx + 1));
    }
}

bool generateAndSaveDiffMap(const std::vector<std::uint8_t> &orig,
                            const std::vector<std::uint8_t> &mod,
                            int width, int height,
                            const std::string &out_path)
{
    size_t px_expected = static_cast<size_t>(width) * static_cast<size_t>(height);
    size_t pxcount = std::min(orig.size(), mod.size());

    if (pxcount == 0) {
        std::cerr << "Warning: zero pixels provided for diff map '" << out_path << "'\n";
        return false;
    }
    if (pxcount != px_expected) {
        std::cerr << "Warning: diff map pixel count (" << pxcount << ") != width*height (" << px_expected << ") for '" << out_path << "'. Using min size.\n";
    }

    int minDiff = std::numeric_limits<int>::max();
    int maxDiff = std::numeric_limits<int>::min();
    std::vector<int> diffs(pxcount);

    for (size_t p = 0; p < pxcount; ++p) {
        int d = static_cast<int>(mod[p]) - static_cast<int>(orig[p]);
        diffs[p] = d;
        if (d < minDiff) minDiff = d;
        if (d > maxDiff) maxDiff = d;
    }

    double s_neg = std::numeric_limits<double>::infinity();
    double s_pos = std::numeric_limits<double>::infinity();

    if (minDiff < 0) {
        s_neg = 127.0 / static_cast<double>(-minDiff);
    }
    if (maxDiff > 0) {
        s_pos = 128.0 / static_cast<double>(maxDiff);
    }

    double s;
    if (!std::isfinite(s_neg) && !std::isfinite(s_pos)) {
        s = 0.0;
    } else {
        s = std::min(s_neg, s_pos);
    }

    std::vector<std::uint8_t> diff_map(px_expected, 127);
    for (size_t p = 0; p < pxcount; ++p) {
        double mapped = 127.0 + s * static_cast<double>(diffs[p]);
        int im = static_cast<int>(std::lround(mapped));
        if (im < 0) im = 0;
        if (im > 255) im = 255;
        diff_map[p] = static_cast<std::uint8_t>(im);
    }

    if (!saveGrayPNG(out_path, diff_map, width, height)) {
        std::cerr << "Failed to save diff map to '" << out_path << "'\n";
        return false;
    }
    return true;
}

int main() {
    const std::string input_path = "images/input.png";

    std::vector<std::uint8_t> secret_pixels;
    int width=0, height=0;
    if (!readGrayPNG(input_path, secret_pixels, width, height)) {
        std::cerr << "Failed to read the input secret image. Exiting.\n";
        return 1;
    }
    size_t secret_size = secret_pixels.size();
    std::cout << "Loaded secret image '" << input_path << "' (" << width << "x" << height << "), bytes: " << secret_size << "\n";

    const unsigned n = 3;
    const unsigned k = 3;
    const unsigned kn = 1;

    const std::string out_base = "results/shamir256test_x";
    fs::create_directories(out_base);
    fs::create_directories(out_base + "/shadows");
    fs::create_directories(out_base + "/shadows_jpeg");
    fs::create_directories(out_base + "/reconstructions");
    fs::create_directories(out_base + "/diff_maps");

    auto field = ss::BinaryField<std::uint8_t>(8);

    // Define the different sets of x_k values to test
    std::vector<std::vector<std::uint8_t>> x_sets = {
        {1, 2, 3},                    // Original baseline
        {1, 255, 254},                // Wraparound set
        {1, 128, 129},                // Boundary adjacent
        {1, 2, 128},                  // Bad result
        {1, 3, 5},                    // Power of primitive
        {85, 170, 255}                // Evenly spaced
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

        // Generate shares with specific x values
        std::vector<std::uint8_t> field_vec = secret_pixels;
        auto shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn, x_values);
        
        if (shares.size() != n) {
            std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
            continue;
        }

        // Save original shares for this set
        for (size_t i = 0; i < shares.size(); ++i) {
            std::string outp = out_base + "/shadows/set" + std::to_string(set_idx+1) + 
                              "_shadow_x" + std::to_string(static_cast<int>(x_values[i])) + ".png";
            saveGrayPNG(outp, shares[i], width, height);
        }

        // Test JPEG compression on each share individually
        for (unsigned share_to_noise = 0; share_to_noise < n; ++share_to_noise) {
            std::cout << "\n--- Noising share " << (share_to_noise + 1) << " (x=" 
                      << static_cast<int>(x_values[share_to_noise]) << ") ---\n";

            // Create a copy of shares where only one is noised
            std::vector<std::vector<std::uint8_t>> comp_shares = shares;

            // Apply JPEG compression to the selected share
            std::string jpeg_path = out_base + "/shadows_jpeg/set" + std::to_string(set_idx+1) + 
                                   "_noised_share" + std::to_string(share_to_noise+1) + 
                                   "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".jpg";

            if (!stbi_write_jpg(jpeg_path.c_str(), width, height, 1, shares[share_to_noise].data(), jpeg_quality)) {
                std::cerr << "Warning: failed to write JPEG for share " << (share_to_noise+1) << " to '" << jpeg_path << "'\n";
                continue;
            }

            // Reload the JPEG compressed share
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

            // Compute metrics for the noised share vs original
            std::cout << "Share " << (share_to_noise+1) << " (x=" << static_cast<int>(x_values[share_to_noise]) 
                      << ") orig vs JPEG | PSNR = " << ss::computePSNR(shares[share_to_noise], comp_shares[share_to_noise]) 
                      << " dB | NPCR = " << ss::computeNPCR(shares[share_to_noise], comp_shares[share_to_noise]) 
                      << "% | UACI = " << ss::computeUACI(shares[share_to_noise], comp_shares[share_to_noise]) << "%\n";

            // Save diff map for the noised share
            std::string diff_path = out_base + "/diff_maps/diff_set" + std::to_string(set_idx+1) + 
                                   "_share" + std::to_string(share_to_noise+1) + 
                                   "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";
            generateAndSaveDiffMap(shares[share_to_noise], comp_shares[share_to_noise], width, height, diff_path);

            // Reconstruct using all three shares (one noised, two clean)
            auto rec = ss::reconstructFromShares<std::uint8_t>(comp_shares, x_values, k, field, kn, static_cast<unsigned>(secret_size));

            // Save reconstruction
            std::string rec_path = out_base + "/reconstructions/reconstruction_set" + std::to_string(set_idx+1) + 
                                  "_noised_share" + std::to_string(share_to_noise+1) + 
                                  "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";
            if (!saveGrayPNG(rec_path, rec, width, height)) {
                std::cerr << "Failed to save reconstruction to '" << rec_path << "'\n";
            }

            // Compute and display reconstruction metrics
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

            // Save diff map for reconstruction
            std::string rec_diff = out_base + "/diff_maps/diff_reconstruction_set" + std::to_string(set_idx+1) + 
                                  "_noised_share" + std::to_string(share_to_noise+1) + 
                                  "_x" + std::to_string(static_cast<int>(x_values[share_to_noise])) + ".png";
            generateAndSaveDiffMap(field_vec, rec, width, height, rec_diff);
        }
    }

    // Save histogram of original secret
    ss::saveHistogram(out_base + "/histogram.dat", ss::computeHistogram(secret_pixels));

    std::cout << "\n=== All tests completed. Total: " << x_sets.size() << " sets × 3 shares = " 
              << (x_sets.size() * 3) << " test cases ===\n";

    return 0;
}