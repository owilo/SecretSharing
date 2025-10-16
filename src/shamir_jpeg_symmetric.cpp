#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <cmath>
#include <limits>

#include "secret.hpp"
#include "image.hpp"

namespace fs = std::filesystem;


static void genCombinations(unsigned center, unsigned r, unsigned start, std::vector<unsigned> &current, std::vector<std::vector<unsigned>> &out) {
    if (current.size() == r) {
        out.push_back(current);
        return;
    }
    for (unsigned d = start; d <= center; ++d) {
        current.push_back(d);
        genCombinations(center, r, d + 1, current, out);
        current.pop_back();
    }
}

static std::string joinIndicesOneBased(const std::vector<unsigned> &indices0based) {
    std::string s;
    for (size_t i = 0; i < indices0based.size(); ++i) {
        if (i) s += "_";
        s += std::to_string(indices0based[i] + 1);
    }
    return s;
}

int main() {
    const std::string input_path = "images/input.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
    size_t secret_size = secret_pixels.size();
    std::cout << "Loaded secret image '" << input_path << "' (" << width << "x" << height << "), bytes: " << secret_size << "\n";

    const unsigned n = 11;
    const unsigned kn = 1;
    const std::string out_base = "results/shamir256jpeg_symmetric";
    fs::create_directories(out_base);
    fs::create_directories(out_base + "/hist");
    fs::create_directories(out_base + "/diff_maps");

    auto field = ss::BinaryField<std::uint8_t>(8);

    std::vector<std::uint8_t> field_vec = secret_pixels;

    ss::saveHistogram(out_base + "/hist/histogram.dat", ss::computeHistogram(field_vec));

    const int jpeg_quality = 100;

    const unsigned center = n / 2;

    const std::vector<unsigned> k_values = {3u, 5u, 7u};

    for (unsigned k : k_values) {
        if ((k % 2) == 0) {
            std::cerr << "Skipping even k=" << k << "\n";
            continue;
        }
        if (k > n) {
            std::cerr << "Skipping k=" << k << " > n=" << n << "\n";
            continue;
        }

        std::string kdir = out_base + "/k" + std::to_string(k);
        fs::create_directories(kdir);
        fs::create_directories(kdir + "/shadows");
        fs::create_directories(kdir + "/shadows_jpeg");
        fs::create_directories(kdir + "/reconstructions");
        fs::create_directories(kdir + "/diff_maps");

        auto shares = ss::getSharesSymmetric<std::uint8_t>(field_vec, k, n, field, kn);
        if (shares.size() != n) {
            std::cerr << "Warning: getSharesSymmetric returned " << shares.size() << " shares (expected " << n << ")\n";
        }

        for (size_t i = 0; i < shares.size(); ++i) {
            std::string outp = kdir + "/shadows/shadow_" + std::to_string(i+1) + ".png";
            ss::saveGrayscalePNG(outp, shares[i], width, height);
        }

        std::vector<std::vector<std::uint8_t>> comp_shares = shares;

        unsigned j = center;
        std::string jpeg_path = kdir + "/shadows_jpeg/shadow_" + std::to_string(j+1) + ".jpg";
        comp_shares[j] = ss::jpegify(shares[j], width, height, jpeg_quality, kdir);

        std::cout << "Share " << (j+1)
            << " orig vs JPEG | PSNR = " << ss::computePSNR(shares[j], comp_shares[j])
            << " dB | NPCR = " << ss::computeNPCR(shares[j], comp_shares[j])
            << "% | UACI = " << ss::computeUACI(shares[j], comp_shares[j])
            << "%\n";

        auto diff_map = ss::generateDiffMap(shares[j], comp_shares[j], width, height);
        std::string diff_path = out_base + "/diff_maps/diff_share_" + std::to_string(j+1) + ".png";
        ss::saveGrayscalePNG(diff_path, diff_map, width, height);

        unsigned r = (k - 1) / 2;
        std::vector<std::vector<unsigned>> combos;
        if (r == 0) {
            combos.push_back(std::vector<unsigned>{});
        } else {
            std::vector<unsigned> cur;
            genCombinations(center, r, 1, cur, combos);
        }

        for (const auto &distances : combos) {
            std::vector<unsigned> indices0;
            indices0.push_back(center);
            for (unsigned d : distances) {
                if (d > center) {
                    std::cerr << "Distance " << d << " out of range for center " << center << "\n";
                    continue;
                }
                indices0.push_back(center - d);
                indices0.push_back(center + d);
            }
            std::sort(indices0.begin(), indices0.end());

            if (indices0.size() != k) {
                std::cerr << "Warning: generated indices size (" << indices0.size() << ") != k (" << k << "). Skipping this combo.\n";
                continue;
            }


            auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices0, comp_shares);

            auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));

            std::string idx_str = joinIndicesOneBased(indices0);
            std::string rec_path = kdir + "/reconstructions/reconstruction_" + idx_str + ".png";
            ss::saveGrayscalePNG(rec_path, rec, width, height);

            std::cout << "k=" << k << " | Reconstruction using shares (" << idx_str << ") | PSNR: " << ss::computePSNR(field_vec, rec) << " dB | NPCR: " << ss::computeNPCR(field_vec, rec) << "% | UACI: " << ss::computeUACI(field_vec, rec) << "%\n";

            std::string rec_diff = kdir + "/diff_maps/diff_reconstruction_" + idx_str + ".png";

            std::vector<std::uint8_t> diff_map = ss::generateDiffMap(field_vec, rec, width, height);
            ss::saveGrayscalePNG(rec_diff, diff_map, width, height);
        }
    }

    return 0;
}
