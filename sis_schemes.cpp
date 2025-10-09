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
    if ((size_t)w * (size_t)h != pixels.size()) {
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

template <typename FieldType>
void runScheme(const std::string &label,
               FieldType field,
               unsigned n,
               unsigned k,
               unsigned kn,
               const std::vector<std::uint8_t> &secret_vec,
               size_t secret_size,
               int width,
               int height,
               const std::string &result_prefix,
               bool do_failed_experiment = false,
               bool divide_share_height_by_k = false)
{
    std::cout << "\n--- " << label << " (n=" << n << ", k=" << k << ", kn=" << kn << ") ---\n";

    fs::create_directories(result_prefix + "/shadows");

    std::vector<std::uint8_t> hist_vec = ss::stretchHistogram(secret_vec, 0, 255, 0, field.getOrder() - 1);

    auto shares = ss::getShares<std::uint8_t>(hist_vec, k, n, field, kn);
    if (shares.size() != n) {
        std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
    }

    int share_h = divide_share_height_by_k ? (height / k) : height;
    for (size_t i = 0; i < shares.size(); ++i) {
        std::string outp = result_prefix + "/shadows/shadow_" + std::to_string(i+1) + ".png";
        if (!saveGrayPNG(outp, shares[i], width, share_h)) {
            std::cerr << "Failed to save " << outp << "\n";
        } else {
            std::cout << "Saved " << outp << "\n";
        }
    }

    std::vector<unsigned> idxA(k);
    std::iota(idxA.begin(), idxA.end(), 0u);

    std::vector<unsigned> idxB;
    for (unsigned i = 2; i < 2 + k && i < shares.size(); ++i) idxB.push_back(i);
    if (idxB.size() < k) {
        for (unsigned i = 0; i < shares.size() && idxB.size() < k; ++i) {
            if (std::find(idxB.begin(), idxB.end(), i) == idxB.end()) idxB.push_back(i);
        }
    }

    auto reconstructAndSave = [&](const std::vector<unsigned> &indices, const std::string &out_path, unsigned rec_k){
        std::vector<std::vector<std::uint8_t>> sel;
        std::vector<std::uint8_t> xs;
        selectSharesAndXs(shares, indices, sel, xs);
        auto rec = ss::reconstructFromShares<std::uint8_t>(sel, xs, rec_k, field, kn, static_cast<unsigned>(secret_size));
        if (!saveGrayPNG(out_path, rec, width, height)) {
            std::cerr << "Failed to save " << out_path << "\n";
        } else {
            std::cout << "Saved " << out_path << "\n";
        }
        return rec;
    };

    std::string recA_path = result_prefix + "/reconstruction1.png";
    auto recA = reconstructAndSave(idxA, recA_path, k);

    std::string recB_path = result_prefix + "/reconstruction2.png";
    auto recB = reconstructAndSave(idxB, recB_path, k);

    std::vector<std::uint8_t> recF;
    std::string recF_path = result_prefix + "/reconstruction_failed.png";
    if (do_failed_experiment) {
        unsigned k_fail = (k > 1 ? k - 1 : 1);
        std::vector<unsigned> idxF;
        for (unsigned i = 0; i < k_fail && i < shares.size(); ++i) idxF.push_back(i);
        recF = reconstructAndSave(idxF, recF_path, k_fail);
    }

    // PSNRs
    std::cout << "\nPSNRs for " << label << " reconstructions:\n";
    try {
        double p1 = ss::computePSNR(secret_vec, recA);
        std::cout << "  " << recA_path << ": " << p1 << " dB\n";
    } catch (...) { std::cout << "  PSNR compute failed for " << recA_path << "\n"; }
    try {
        double p2 = ss::computePSNR(secret_vec, recB);
        std::cout << "  " << recB_path << ": " << p2 << " dB\n";
    } catch (...) { std::cout << "  PSNR compute failed for " << recB_path << "\n"; }
    if (do_failed_experiment) {
        try {
            double pf = ss::computePSNR(secret_vec, recF);
            std::cout << "  " << recF_path << ": " << pf << " dB\n";
        } catch (...) { std::cout << "  PSNR compute failed for " << recF_path << "\n"; }
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

    // 1) Shamir GF(251) (kn = 1)
    runScheme("Shamir GF(251)", ss::PrimeField<std::uint8_t>(251), n, k, 1, secret_pixels, secret_size, width, height, "results/shamir251", true, false);

    // 2) Shamir GF(2^8) (AES poly 0x11B) kn = 1
    runScheme("Shamir GF(2^8)", ss::BinaryField<std::uint8_t>(8), n, k, 1, secret_pixels, secret_size, width, height, "results/shamir256", true, false);

    // 3) Thien-Lin GF(251) (kn = 4)
    runScheme("Thien-Lin GF(251)", ss::PrimeField<std::uint8_t>(251), n, k, 4, secret_pixels, secret_size, width, height, "results/thienlin251", false, true);

    // 4) Thien-Lin GF(2^8) (kn = 4)
    runScheme("Thien-Lin GF(2^8)", ss::BinaryField<std::uint8_t>(8), n, k, 4, secret_pixels, secret_size, width, height, "results/thienlin256", false, true);
    return 0;
}