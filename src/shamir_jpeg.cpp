#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <limits>

#include "ss/secret.hpp"
#include "ss/image.hpp"

namespace fs = std::filesystem;

int main() {
    const std::string input_path = "images/baboon.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
    size_t secret_size = secret_pixels.size();
    std::cout << "Loaded secret image '" << input_path << "' (" << width << "x" << height << "), bytes: " << secret_size << "\n";

    const unsigned n = 5;
    const unsigned k = 3;
    const unsigned kn = 1;

    const std::string out_base = "results/shamir256jpeg";
    fs::create_directories(out_base);
    fs::create_directories(out_base + "/shadows");
    fs::create_directories(out_base + "/shadows_jpeg");
    fs::create_directories(out_base + "/reconstructions");
    fs::create_directories(out_base + "/diff_maps");

    auto field = ss::BinaryField<std::uint8_t>(8);

    // may be useless here
    //std::vector<std::uint8_t> field_vec = ss::clampPixels(secret_pixels, 0, static_cast<std::uint8_t>(field.getOrder()-1));
    std::vector<std::uint8_t> field_vec = secret_pixels;

    auto shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn);
    if (shares.size() != n) {
        std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
    }

    for (size_t i = 0; i < shares.size(); ++i) {
        std::string outp = out_base + "/shadows/shadow_" + std::to_string(i+1) + ".png";
        ss::saveGrayscalePNG(outp, shares[i], width, height);
    }

    ss::saveHistogram(out_base + "/histogram.dat", ss::computeHistogram(field_vec));

    const int jpeg_quality = 100;

    std::vector<std::vector<std::uint8_t>> comp_shares = shares;

    std::vector<unsigned> noisy_indices = {2u, 3u, 4u};

    for (unsigned j_idx = 0; j_idx < noisy_indices.size(); ++j_idx) {
        unsigned j = noisy_indices[j_idx];

        std::string jpeg_path = out_base + "/shadows_jpeg/shadow_" + std::to_string(j+1) + ".jpg";
        comp_shares[j] = ss::jpegify(shares[j], width, height, jpeg_quality, jpeg_path);

        std::cout << "Share " << (j+1)
                << " orig vs JPEG | PSNR = " << ss::computePSNR(shares[j], comp_shares[j])
                << " dB | NPCR = " << ss::computeNPCR(shares[j], comp_shares[j])
                << "% | UACI = " << ss::computeUACI(shares[j], comp_shares[j])
                << "%\n";

        auto diff_map = ss::generateDiffMap(shares[j], comp_shares[j], width, height);
        std::string diff_path = out_base + "/diff_maps/diff_share_" + std::to_string(j+1) + ".png";
        ss::saveGrayscalePNG(diff_path, diff_map, width, height);
    }

    std::vector<std::array<unsigned,3>> combos = {
        {0u, 1u, 2u},
        {0u, 1u, 3u},
        {0u, 1u, 4u},
        {0u, 3u, 4u},
        {1u, 3u, 4u}
    };

    for (const auto &combo : combos) {
        std::vector<unsigned> indices = {combo[0], combo[1], combo[2]};
        auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, comp_shares);

        auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));

        std::string rec_path = out_base + "/reconstructions/reconstruction_";
        rec_path += std::to_string(indices[0]+1) + "_" + std::to_string(indices[1]+1) + "_" + std::to_string(indices[2]+1) + ".png";
        ss::saveGrayscalePNG(rec_path, rec, width, height);

        std::cout << "Reconstruction using shares (" << (indices[0]+1) << "," << (indices[1]+1) << "," << (indices[2]+1) << ") | PSNR: " << ss::computePSNR(field_vec, rec) << " dB | NPCR: " << ss::computeNPCR(field_vec, rec) << "% | UACI: " << ss::computeUACI(field_vec, rec) << "%\n";

        std::string rec_diff = out_base + "/diff_maps/diff_reconstruction_" + std::to_string(indices[0]+1) + "_" + std::to_string(indices[1]+1) + "_" + std::to_string(indices[2]+1) + ".png";
        std::vector<std::uint8_t> diff_map = ss::generateDiffMap(field_vec, rec, width, height);
        ss::saveGrayscalePNG(rec_diff, diff_map, width, height);
    }

    return 0;
}