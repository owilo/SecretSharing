#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <numeric>

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

    const unsigned n = 6;
    const unsigned k = 4;
    const unsigned kn = 1;

    const std::string out_base = "results/shamir256jpeg";
    fs::create_directories(out_base);
    fs::create_directories(out_base + "/shadows");
    fs::create_directories(out_base + "/shadows_jpeg");
    fs::create_directories(out_base + "/reconstructions");

    auto field = ss::BinaryField<std::uint8_t>(8);

    std::vector<std::uint8_t> field_vec = ss::clampPixels(secret_pixels, 0, static_cast<std::uint8_t>(field.getOrder()-1));

    auto shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn);
    if (shares.size() != n) {
        std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
    }

    for (size_t i = 0; i < shares.size(); ++i) {
        std::string outp = out_base + "/shadows/shadow_" + std::to_string(i+1) + ".png";
        saveGrayPNG(outp, shares[i], width, height);
    }

    ss::saveHistogram(out_base + "/histogram.dat", ss::computeHistogram(field_vec));

    const int jpeg_quality = 98;

    for (unsigned j = 0; j < n; ++j) {
        std::vector<std::vector<std::uint8_t>> comp_shares = shares;

        std::string jpeg_path = out_base + "/shadows_jpeg/shadow_" + std::to_string(j+1) + ".jpg";
        if (!stbi_write_jpg(jpeg_path.c_str(), width, height, 1, comp_shares[j].data(), jpeg_quality)) {
            std::cerr << "Warning: failed to write JPEG for share " << j+1 << " to '" << jpeg_path << "'\n";
        } else {
            int w2=0, h2=0, channels=0;
            unsigned char *data = stbi_load(jpeg_path.c_str(), &w2, &h2, &channels, 1);
            if (!data) {
                std::cerr << "Warning: failed to reload JPEG '" << jpeg_path << "': " << stbi_failure_reason() << "\n";
            } else {
                if (w2 != width || h2 != height) {
                    std::cerr << "Warning: reloaded JPEG size mismatch for '" << jpeg_path << "' (" << w2 << "x" << h2 << ") expected (" << width << "x" << height << ")\n";
                }
                comp_shares[j].assign(data, data + (w2 * h2));
                stbi_image_free(data);
            }
        }

        std::cout << "Share " << (j+1) << " orig vs JPEG | "
                  << "PSNR = " << ss::computePSNR(shares[j], comp_shares[j]) << " dB"
                  << " | NPCR = " << ss::computeNPCR(shares[j], comp_shares[j]) << "%"
                  << " | UACI = " << ss::computeUACI(shares[j], comp_shares[j]) << "%\n";

        std::vector<unsigned> indices;
        for (unsigned t = 0; t < k; ++t) {
            indices.push_back(static_cast<unsigned>((j + t) % n));
        }

        std::vector<std::vector<std::uint8_t>> sel;
        std::vector<std::uint8_t> xs;
        selectSharesAndXs(comp_shares, indices, sel, xs);

        auto rec = ss::reconstructFromShares<std::uint8_t>(sel, xs, k, field, kn, static_cast<unsigned>(secret_size));

        std::string rec_path = out_base + "/reconstructions/reconstruction_share_" + std::to_string(j+1) + ".png";
        if (!saveGrayPNG(rec_path, rec, width, height)) {
            std::cerr << "Failed to save reconstruction to '" << rec_path << "'\n";
        }

        std::cout << "Recon share " << (j+1) << " | "
                  << " | PSNR: " << ss::computePSNR(field_vec, rec) << " dB"
                  << " | NPCR: " << ss::computeNPCR(field_vec, rec) << "%"
                  << " | UACI: " << ss::computeUACI(field_vec, rec) << "%\n";
    }

    std::cout << "All reconstructions saved under: " << out_base << "/reconstructions/ (raw shadows in " << out_base << "/shadows and compressed share files in " << out_base << "/shadows_jpeg)\n";

    return 0;
}
