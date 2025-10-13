#include "image.hpp"

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <fstream>

namespace ss {

std::vector<std::uint8_t> stretchHistogram(
    const std::vector<std::uint8_t>& image,
    std::uint8_t in_min,
    std::uint8_t in_max,
    std::uint8_t out_min,
    std::uint8_t out_max
) noexcept {
    const int inMin  = static_cast<int>(in_min);
    const int inMax  = static_cast<int>(in_max);
    const int outMin = static_cast<int>(out_min);
    const int outMax = static_cast<int>(out_max);

    if (inMin == outMin && inMax == outMax) {
        return image;
    }

    if (inMax == inMin) {
        return std::vector<std::uint8_t>(image.size(), outMax);
    }

    const double scale = static_cast<double>(outMax - outMin) / static_cast<double>(inMax - inMin);

    std::vector<std::uint8_t> result;
    result.reserve(image.size());

    for (const auto &p : image) {
        const int pv = static_cast<int>(p);
        std::uint8_t newValue;
        if (pv <= inMin) {
            newValue = static_cast<std::uint8_t>(std::clamp(outMin, 0, 255));
        } else if (pv >= inMax) {
            newValue = static_cast<std::uint8_t>(std::clamp(outMax, 0, 255));
        } else {
            const double mapped = outMin + (pv - inMin) * scale;
            const int rounded = static_cast<int>(std::lround(mapped));
            newValue = static_cast<std::uint8_t>(std::clamp(rounded, 0, 255));
        }
        result.push_back(newValue);
    }

    return result;
}

std::array<unsigned, 256> computeHistogram(const std::vector<std::uint8_t>& image) {
    std::array<unsigned, 256> hist = {};
    for (const auto &p : image) {
        ++hist[p];
    }
    return hist;
}

bool saveHistogram(std::string path, const std::array<unsigned, 256>& histogram) {
    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        return false;
    }

    for (size_t i = 0; i < histogram.size(); ++i) {
        ofs << i << " " << histogram[i] << "\n";
    }
    return true;
}

double computeEntropyPerPixel(const std::vector<std::uint8_t>& image) {
    auto hist = computeHistogram(image);
    const double total = static_cast<double>(image.size());
    if (total == 0.0) {
        return 0.0;
    }

    double entropy = 0.0;
    for (const auto &count : hist) {
        if (count > 0) {
            double p = static_cast<double>(count) / total;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

std::vector<std::uint8_t> clampPixels(const std::vector<std::uint8_t>& image, std::uint8_t min_val, std::uint8_t max_val) {
    std::vector<std::uint8_t> result;
    result.reserve(image.size());
    for (const auto &p : image) {
        result.push_back(std::clamp(p, min_val, max_val));
    }
    return result;
}

double computePSNR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon) {
    if (orig.size() != recon.size()) {
        throw std::invalid_argument("Vectors must be the same size for PSNR computation");
    }
    double mse = 0.0;
    const std::size_t n = orig.size();
    for (std::size_t i = 0; i < n; ++i) {
        double diff = static_cast<double>(orig[i]) - static_cast<double>(recon[i]);
        mse += diff * diff;
    }
    mse /= n;
    if (mse == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return 10.0 * std::log10(65025.0 / mse);
}

double computeNPCR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon) {
    if (orig.size() != recon.size()) {
        throw std::invalid_argument("Vectors must be the same size for NPCR computation");
    }
    std::size_t diff_count = 0;
    const std::size_t n = orig.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (orig[i] != recon[i]) {
            ++diff_count;
        }
    }
    return (static_cast<double>(diff_count) / static_cast<double>(n)) * 100.0;
}

double computeUACI(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon) {
    if (orig.size() != recon.size()) {
        throw std::invalid_argument("Vectors must be the same size for UACI computation");
    }
    double total_diff = 0.0;
    const std::size_t n = orig.size();
    for (std::size_t i = 0; i < n; ++i) {
        total_diff += std::abs(static_cast<double>(orig[i]) - static_cast<double>(recon[i]));
    }
    return (total_diff / (static_cast<double>(n) * 255.0)) * 100.0;
}

} // namespace ss