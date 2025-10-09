#pragma once

#include <array>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace ss {

std::vector<std::uint8_t> stretchHistogram(const std::vector<std::uint8_t>& image, std::uint8_t in_min, std::uint8_t in_max, std::uint8_t out_min, std::uint8_t out_max) noexcept;

std::array<unsigned, 256> computeHistogram(const std::vector<std::uint8_t>& image);

std::vector<std::uint8_t> clampPixels(const std::vector<std::uint8_t>& image, std::uint8_t min_val, std::uint8_t max_val);

bool saveHistogram(std::string path, const std::array<unsigned, 256>& histogram);

double computePSNR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeNPCR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeUACI(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

} // namespace ss