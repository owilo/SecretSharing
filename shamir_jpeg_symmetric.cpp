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


bool generateAndSaveDiffMap(const std::vector<std::uint8_t> &orig, const std::vector<std::uint8_t> &mod, int width, int height, const std::string &out_path) {
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

        {
            unsigned j = center;
            std::string jpeg_path = kdir + "/shadows_jpeg/shadow_" + std::to_string(j+1) + ".jpg";

            if (!stbi_write_jpg(jpeg_path.c_str(), width, height, 1, shares[j].data(), jpeg_quality)) {
                std::cerr << "Warning: failed to write JPEG for share " << (j+1) << " to '" << jpeg_path << "'\n";
            } else {
                int w2=0, h2=0, channels=0;
                unsigned char *data = stbi_load(jpeg_path.c_str(), &w2, &h2, &channels, 1);
                if (!data) {
                    std::cerr << "Warning: failed to reload JPEG '" << jpeg_path << "': " << stbi_failure_reason() << "\n";
                } else {
                    if (w2 != width || h2 != height) {
                        std::cerr << "Warning: reloaded JPEG size mismatch for '" << jpeg_path << "' (" << w2 << "x" << h2 << ") expected (" << width << "x" << height << ")\n";
                    }
                    comp_shares[j].assign(data, data + (static_cast<size_t>(w2) * static_cast<size_t>(h2)));
                    stbi_image_free(data);

                    std::cout << "k=" << k << " | Share " << (j+1) << " orig vs JPEG | PSNR = " << ss::computePSNR(shares[j], comp_shares[j]) << " dB | NPCR = " << ss::computeNPCR(shares[j], comp_shares[j]) << "% | UACI = " << ss::computeUACI(shares[j], comp_shares[j]) << "%\n";

                    std::string diff_path = kdir + "/diff_maps/diff_share_" + std::to_string(j+1) + ".png";
                    generateAndSaveDiffMap(shares[j], comp_shares[j], width, height, diff_path);
                }
            }
        }

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
            generateAndSaveDiffMap(field_vec, rec, width, height, rec_diff);
        }
    }

    return 0;
}
