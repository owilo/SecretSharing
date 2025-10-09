#pragma once

#include <vector>
#include <cstdint>
#include <stdexcept>

namespace ss {

std::vector<std::uint8_t> stretchHistogram(const std::vector<std::uint8_t>& image, std::uint8_t in_min, std::uint8_t in_max, std::uint8_t out_min, std::uint8_t out_max) noexcept;

double computePSNR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeNPCR(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

double computeUACI(const std::vector<std::uint8_t>& orig, const std::vector<std::uint8_t>& recon);

} // namespace ss