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

    ss::saveGrayscalePNG(out_path, diff_map, width, height);
    return true;
}

int main() {
    const std::string input_path = "images/input.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
    size_t secret_size = secret_pixels.size();
    std::cout << "Loaded secret image '" << input_path << "' (" << width << "x" << height << "), bytes: " << secret_size << "\n";

    const unsigned n = 6;
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

    std::vector<unsigned> noisy_indices = {2u, 3u, 4u, 5u};

    for (unsigned j_idx = 0; j_idx < noisy_indices.size(); ++j_idx) {
        unsigned j = noisy_indices[j_idx];
        std::string jpeg_path = out_base + "/shadows_jpeg/shadow_" + std::to_string(j+1) + ".jpg";

        if (!stbi_write_jpg(jpeg_path.c_str(), width, height, 1, shares[j].data(), jpeg_quality)) {
            std::cerr << "Warning: failed to write JPEG for share " << (j+1) << " to '" << jpeg_path << "'\n";
            continue;
        }

        int w2=0, h2=0, channels=0;
        unsigned char *data = stbi_load(jpeg_path.c_str(), &w2, &h2, &channels, 1);
        if (!data) {
            std::cerr << "Warning: failed to reload JPEG '" << jpeg_path << "': " << stbi_failure_reason() << "\n";
            continue;
        } else {
            if (w2 != width || h2 != height) {
                std::cerr << "Warning: reloaded JPEG size mismatch for '" << jpeg_path << "' (" << w2 << "x" << h2 << ") expected (" << width << "x" << height << ")\n";
            }
            comp_shares[j].assign(data, data + (static_cast<size_t>(w2) * static_cast<size_t>(h2)));
            stbi_image_free(data);

            std::cout << "Share " << (j+1) << " orig vs JPEG | PSNR = " << ss::computePSNR(shares[j], comp_shares[j]) << " dB | NPCR = " << ss::computeNPCR(shares[j], comp_shares[j]) << "% | UACI = " << ss::computeUACI(shares[j], comp_shares[j]) << "%\n";

            std::string diff_path = out_base + "/diff_maps/diff_share_" + std::to_string(j+1) + ".png";
            generateAndSaveDiffMap(shares[j], comp_shares[j], width, height, diff_path);
        }
    }

    for (unsigned noisy : noisy_indices) {
        std::vector<unsigned> indices = {0u, 1u, noisy};
        auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, comp_shares);

        auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));

        std::string rec_path = out_base + "/reconstructions/reconstruction_0_1_" + std::to_string(noisy+1) + ".png";
        ss::saveGrayscalePNG(rec_path, rec, width, height);

        std::cout << "Reconstruction using shares (0,1," << (noisy+1) << ") | PSNR: " << ss::computePSNR(field_vec, rec) << " dB | NPCR: " << ss::computeNPCR(field_vec, rec) << "% | UACI: " << ss::computeUACI(field_vec, rec) << "%\n";

        std::string rec_diff = out_base + "/diff_maps/diff_reconstruction_0_1_" + std::to_string(noisy+1) + ".png";
        generateAndSaveDiffMap(field_vec, rec, width, height, rec_diff);
    }

    return 0;
}
