#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <filesystem>

#include "secret.hpp"
#include "image.hpp"

namespace fs = std::filesystem;

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
               bool divide_share_height_by_k = false,
               bool clampToField = true)
{
    std::cout << "\n--- " << label << " (n=" << n << ", k=" << k << ", kn=" << kn << ") ---\n";

    fs::create_directories(result_prefix + "/shadows");

    std::vector<std::uint8_t> field_vec;
    if (clampToField) {
        field_vec = ss::clampPixels(secret_vec, 0, static_cast<std::uint8_t>(field.getOrder() - 1));
    } else {
        field_vec = ss::stretchHistogram(secret_vec, 0, 255, 0, field.getOrder() - 1);
    }

    auto shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn);
    if (shares.size() != n) {
        std::cerr << "Warning: getShares returned " << shares.size() << " shares (expected " << n << ")\n";
    }

    int share_h = divide_share_height_by_k ? (height / k) : height;
    for (size_t i = 0; i < shares.size(); ++i) {
        std::string outp = result_prefix + "/shadows/shadow_" + std::to_string(i+1) + ".png";
        ss::saveGrayscalePNG(outp, shares[i], width, share_h);
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

    auto reconstructAndSave = [&](const std::vector<unsigned>& indices, const std::string& out_path, unsigned rec_k){
        auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, shares);
        auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, rec_k, field, kn, static_cast<unsigned>(secret_size));
        ss::saveGrayscalePNG(out_path, rec, width, height);
        return rec;
    };

    ss::saveHistogram(result_prefix + "/histogram.dat", ss::computeHistogram(field_vec));

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

    std::cout << "  " << recA_path << " | PSNR : " << ss::computePSNR(secret_vec, recA) << " dB | NPCR : " << ss::computeNPCR(secret_vec, recA) << "% | UACI : " << ss::computeUACI(secret_vec, recA) << "%\n";
    std::cout << "  " << recB_path << " | PSNR : " << ss::computePSNR(secret_vec, recB) << " dB | NPCR : " << ss::computeNPCR(secret_vec, recB) << "% | UACI : " << ss::computeUACI(secret_vec, recB) << "%\n";
    if (do_failed_experiment) {
        std::cout << "  " << recF_path << " | PSNR : " << ss::computePSNR(secret_vec, recF) << " dB | NPCR : " << ss::computeNPCR(secret_vec, recF) << "% | UACI : " << ss::computeUACI(secret_vec, recF) << "%\n";
    }

    if (kn == 1) {
        std::string share1_path = result_prefix + "/shadows/shadow_1.png";
        const auto &share1 = shares[0];
        std::cout << "  " << share1_path << " | PSNR : " << ss::computePSNR(secret_vec, share1) << " dB | NPCR: " << ss::computeNPCR(secret_vec, share1) << "% | UACI: " << ss::computeUACI(secret_vec, share1) << "%\n";
    }

}

int main() {
    const std::string input_path = "images/input.png";
    auto [secret_pixels, width, height] = ss::readGrayscalePNG(input_path);
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