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

void printMetrics(const std::vector<std::uint8_t>& image1, const std::vector<std::uint8_t>& image2, std::string image1name, std::string image2name) {
    std::cout << image1name << " vs " << image2name 
        << " : orig vs JPEG | PSNR = " << ss::computePSNR(image1, image2)
        << " dB | NPCR = " << ss::computeNPCR(image1, image2)
        << "% | UACI = " << ss::computeUACI(image1, image2)
        << "%\n";
}

int main() {
    const std::string input_path = "images/input.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
    size_t secret_size = secret_pixels.size();

    const unsigned n = 11;
    const unsigned k = 3;
    const unsigned kn = 1;
    std::vector<unsigned> indices = {3, 6, 9};
    const int jpeg_quality = 100;

    const std::string out_base = "results/shamir256jpeg_median";
    fs::create_directories(out_base);

    auto field = ss::BinaryField<std::uint8_t>(8);

    std::vector<std::uint8_t> input = secret_pixels;

    auto shares = ss::getShares<std::uint8_t>(input, k, n, field, kn);

    auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, shares);

    auto& sj0 = selectedShares[1];
    std::string jpeg_path = out_base + "/share_jpeg.jpg";
    sj0 = ss::jpegify(sj0, width, height, jpeg_quality, jpeg_path);

    std::vector<std::uint8_t> jsm1(sj0.size());
    std::transform(sj0.begin(), sj0.end(), jsm1.begin(), [](int x) { return x == 0 ? x : x - 1; });
    ss::saveGrayscalePNG(out_base + "/share_jpeg_minus1.png", jsm1, width, height);

    std::vector<std::uint8_t> jsp1(sj0.size());
    std::transform(sj0.begin(), sj0.end(), jsp1.begin(), [](int x) { return x == 255 ? x : x + 1; });
    ss::saveGrayscalePNG(out_base + "/share_jpeg_plus1.png", jsp1, width, height);

    auto ir0 = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));
    ss::saveGrayscalePNG(out_base + "/image_jpeg.png", ir0, width, height);

    const std::vector<std::vector<std::uint8_t>> shares_m1 = {selectedShares[0], jsm1, selectedShares[2]};
    const std::vector<std::vector<std::uint8_t>> shares_p1 = {selectedShares[0], jsp1, selectedShares[2]};

    auto irm1 = ss::reconstructFromShares<std::uint8_t>(shares_m1, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));
    ss::saveGrayscalePNG(out_base + "/image_jpeg_minus1.png", irm1, width, height);

    auto irp1 = ss::reconstructFromShares<std::uint8_t>(shares_p1, selectedXs, k, field, kn, static_cast<unsigned>(secret_size));
    ss::saveGrayscalePNG(out_base + "/image_jpeg_plus1.png", irp1, width, height);

    auto if0 = ss::medianFilter(ir0, width, height, 3);
    ss::saveGrayscalePNG(out_base + "/image_median.png", if0, width, height);

    std::vector<std::uint8_t> ic(ir0.size());
    for (std::size_t i = 0; i < ic.size(); ++i) {
        std::uint8_t v = if0[i];
        std::uint8_t candidates[3] = { ir0[i], irm1[i], irp1[i] };
        std::uint8_t best = candidates[0];
        int min_diff = std::abs(int(candidates[0]) - int(v));
        for (int j = 1; j < 3; ++j) {
            int diff = std::abs(int(candidates[j]) - int(v));
            if (diff < min_diff) {
                min_diff = diff;
                best = candidates[j];
            }
        }
        ic[i] = best;
    }
    ss::saveGrayscalePNG(out_base + "/image_corrected.png", ic, width, height);

    printMetrics(input, if0, "original", "image_median");
    printMetrics(input, ir0, "original", "image_jpeg");
    printMetrics(input, irm1, "original", "image_jpeg_minus1");
    printMetrics(input, irp1, "original", "image_jpeg_plus1");
    printMetrics(input, ic, "original", "image_corrected");

    std::vector<std::uint8_t> diff_map_if0 = ss::generateDiffMap(input, if0, width, height);
    ss::saveGrayscalePNG(out_base + "/diff_original_median.png", diff_map_if0, width, height);

    std::vector<std::uint8_t> diff_map_ic = ss::generateDiffMap(input, ic, width, height);
    ss::saveGrayscalePNG(out_base + "/diff_original_corrected.png", diff_map_ic, width, height);

    return 0;
}
