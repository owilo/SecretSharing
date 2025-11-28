#pragma once

#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <string_view>
#include <array>
#include <tuple>
#include <stdexcept>
#include <sstream>
#include <span>

extern "C" {
    #include "stb_image.h"
    #include "stb_image_write.h"
}

namespace ss {

std::pair<int, int> saveGrayscalePNG(const std::filesystem::path& path, std::span<const std::uint8_t> pixels, int w, int h);

std::tuple<std::vector<std::uint8_t>, int, int> readGrayscale(const std::filesystem::path& path);

std::vector<std::uint8_t> generateDiffMap(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon, int width, int height);
    
std::vector<std::uint8_t> stretchHistogram(const std::vector<std::uint8_t>& image, std::uint8_t in_min, std::uint8_t in_max, std::uint8_t out_min, std::uint8_t out_max) noexcept;

std::array<unsigned, 256> computeHistogram(const std::vector<std::uint8_t>& image);

std::vector<std::uint8_t> medianFilter(const std::vector<std::uint8_t>& image, int width, int height, int filter_size, bool circular = false);

std::vector<std::uint8_t> jpegify(const std::vector<std::uint8_t>& image, int width, int height, int quality, std::string output_file = "");

double computeEntropyPerPixel(const std::vector<std::uint8_t>& image);

std::vector<std::uint8_t> clampPixels(const std::vector<std::uint8_t>& image, std::uint8_t min_val, std::uint8_t max_val);

bool saveHistogram(std::string path, const std::array<unsigned, 256>& histogram);

double computePSNR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeNPCR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeUACI(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeEntropy(const std::vector<std::uint8_t>& image);

} // namespace ss