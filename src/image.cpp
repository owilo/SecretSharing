#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "image.hpp"

namespace ss {

std::pair<int, int> saveGrayscalePNG(const std::filesystem::path& path, std::span<const std::uint8_t> pixels, int w, int h) {
    if (w <= 0 || h <= 0) {
        throw std::invalid_argument("saveGrayPNG: width and height must be positive");
    }

    const std::size_t expected = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
    if (pixels.size() != expected) {
        std::ostringstream oss;
        oss << "saveGrayPNG: pixel buffer size (" << pixels.size() << ") != w * h (" << expected << ")";
        throw std::invalid_argument(oss.str());
    }

    if (stbi_write_png(path.string().c_str(), w, h, 1, pixels.data(), w) == 0) {
        throw std::runtime_error("saveGrayPNG: stbi_write_png failed for " + path.string());
    }

    return {w, h};
}

std::tuple<std::vector<std::uint8_t>, int, int> readGrayscalePNG(const std::filesystem::path& path) {
    int w = 0, h = 0, channels = 0;
    std::uint8_t *data = stbi_load(path.string().c_str(), &w, &h, &channels, 1);
    if (!data) {
        std::string reason = stbi_failure_reason() ? stbi_failure_reason() : "unknown";
        throw std::runtime_error("readGrayPNG: failed to load '" + path.string() + "': " + reason);
    }

    const std::size_t expected = static_cast<std::size_t>(w) * static_cast<std::size_t>(h);
    std::vector<std::uint8_t> pixels;
    pixels.assign(data, data + expected);

    stbi_image_free(data);
    return {std::move(pixels), w, h};
}

std::vector<std::uint8_t> generateDiffMap(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon, int width, int height) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("generateDiffMap: width and height must be positive");
    }

    const std::size_t px_expected = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);

    if (orig.size() != px_expected || recon.size() != px_expected) {
        throw std::invalid_argument("generateDiffMap: input buffers must have exactly width*height elements");
    }

    int minDiff = std::numeric_limits<int>::max();
    int maxDiff = std::numeric_limits<int>::min();

    std::vector<int> diffs;
    diffs.reserve(px_expected);

    for (std::size_t i = 0; i < px_expected; ++i) {
        int d = static_cast<int>(recon[i]) - static_cast<int>(orig[i]);
        diffs.push_back(d);
        if (d < minDiff) minDiff = d;
        if (d > maxDiff) maxDiff = d;
    }

    // Determine scaling
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

    std::vector<std::uint8_t> diff_map(px_expected, 127u);

    for (std::size_t i = 0; i < px_expected; ++i) {
        double mapped = 127.0 + s * static_cast<double>(diffs[i]);
        int im = static_cast<int>(std::lround(mapped));
        if (im < 0) im = 0;
        else if (im > 255) im = 255;
        diff_map[i] = static_cast<std::uint8_t>(im);
    }

    return diff_map;
}

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

std::vector<std::uint8_t> medianFilter(const std::vector<std::uint8_t>& image, int width, int height, int filter_size, bool circular) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("Width and height must be positive");
    }
    if (filter_size <= 0 || filter_size % 2 == 0) {
        throw std::invalid_argument("Filter_size must be a positive odd integer");
    }

    const std::size_t px_expected = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (image.size() != px_expected) {
        throw std::invalid_argument("Input buffer must have exactly width*height elements");
    }

    const int radius = filter_size / 2;
    std::vector<std::uint8_t> filtered(image.size(), 0);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            std::vector<std::uint8_t> neighborhood;
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    if (circular && (dx * dx + dy * dy > radius * radius)) {
                        continue;
                    }
                    int nx = x + dx;
                    int ny = y + dy;
                    if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                        neighborhood.push_back(image[ny * width + nx]);
                    }
                }
            }
            if (!neighborhood.empty()) {
                std::nth_element(neighborhood.begin(), neighborhood.begin() + neighborhood.size() / 2, neighborhood.end());
                filtered[y * width + x] = neighborhood[neighborhood.size() / 2];
            } else {
                filtered[y * width + x] = image[y * width + x];
            }
        }
    }

    return filtered;
}

// TODO : use in-memory operations if possible
std::vector<std::uint8_t> jpegify(const std::vector<std::uint8_t>& image, int width, int height, int quality, std::string output_file) {
    if (image.empty()) {
        throw std::invalid_argument("jInput image is empty");
    }
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("Width and height must be positive");
    }
    if (quality < 1 || quality > 100) {
        throw std::invalid_argument("Quality must be in [1,100]");
    }
    const std::size_t expected = static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
    if (image.size() != expected) {
        throw std::invalid_argument("Input image size does not match width*height");
    }

    const bool use_temp = output_file.empty();
    std::string path = use_temp ? "tmp/tmp.jpg" : output_file;

    if (!stbi_write_jpg(path.c_str(), width, height, 1, image.data(), quality)) {
        if (use_temp) { std::error_code ec; std::filesystem::remove(path, ec); }
        throw std::runtime_error("Failed to write JPEG to '" + path + "'");
    }

    int w2 = 0, h2 = 0, channels = 0;
    unsigned char* reloaded = stbi_load(path.c_str(), &w2, &h2, &channels, 1);
    if (!reloaded) {
        if (use_temp) { std::error_code ec; std::filesystem::remove(path, ec); }
        throw std::runtime_error(std::string("Failed to reload JPEG: ") + stbi_failure_reason());
    }

    if (w2 != width || h2 != height) {
        stbi_image_free(reloaded);
        if (use_temp) { std::error_code ec; std::filesystem::remove(path, ec); }
        throw std::runtime_error("jpegify: reloaded JPEG size mismatch");
    }

    std::vector<std::uint8_t> result;
    result.reserve(static_cast<std::size_t>(w2) * static_cast<std::size_t>(h2));
    result.insert(result.end(), reloaded, reloaded + (static_cast<std::size_t>(w2) * static_cast<std::size_t>(h2)));

    stbi_image_free(reloaded);

    if (use_temp) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    return result;
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