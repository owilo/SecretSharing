#include <vector>
#include <cstdint>
#include <stdexcept>
#include <random>
#include <iostream>
#include <fstream>
#include <cmath>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

using data_t = std::vector<std::uint8_t>;

double computePSNR(const std::vector<uint8_t>& orig, const std::vector<uint8_t>& recon) {
    if (orig.size() != recon.size()) {
        throw std::invalid_argument("Vectors must be the same size for PSNR computation");
    }

    double mse = 0.0;
    size_t N = orig.size();
    for (size_t i = 0; i < N; ++i) {
        int diff = static_cast<int>(orig[i]) - static_cast<int>(recon[i]);
        mse += diff * diff;
    }

    mse /= static_cast<double>(N);
    if (mse == 0) return INFINITY;

    double psnr = 10.0 * log10((255.0 * 255.0) / mse);
    return psnr;
}

std::uint8_t evaluatePolynomial(const data_t& coefficients, std::uint8_t x, unsigned field) {
    int result = 0;
    int xPower = 1; // x^0
    int xi = static_cast<int>(x);
    int m = static_cast<int>(field);

    for (const auto& coeff_u8 : coefficients) {
        int coeff = static_cast<int>(coeff_u8);
        result = (result + coeff * xPower) % m;
        xPower = (xPower * xi) % m;
    }

    if (result < 0) result += m;
    return static_cast<std::uint8_t>(result);
}

std::vector<data_t> getShares(const data_t& data, unsigned k, unsigned n, unsigned field, unsigned kn = 1) {
    // Allocate space for shares
    unsigned shareSize = (data.size() + kn - 1) / kn;

    std::vector<data_t> shares(n);
    for (unsigned i = 0; i < n; ++i) {
        shares[i].reserve(shareSize);
    }

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> random(0, field - 1);

    unsigned dataCursor = 0;
    for (unsigned i = 0; i < shareSize; ++i) {
        std::vector<std::uint8_t> polynomialCoefficients(k);
        // Fill kn bytes of data into the polynomial coefficients
        for (unsigned j = 0; j < kn; ++j) {
            if (dataCursor < data.size()) {
                polynomialCoefficients[j] = std::min(data[dataCursor++], static_cast<std::uint8_t>(field - 1));
            } else {
                polynomialCoefficients[j] = 0;
            }
        }
        // Fill the rest with random coefficients
        for (unsigned j = kn; j < k; ++j) {
            polynomialCoefficients[j] = random(rng) % field;
        }

        // Calculate shares
        for (unsigned j = 0; j < n; ++j) {
            shares[j].push_back(evaluatePolynomial(polynomialCoefficients, j + 1, field));
        }
    }

    return shares;
}

static int mod(int a, int m) {
    int r = a % m;
    if (r < 0) r += m;
    return r;
}

static int modInverseInt(int a, int m) {
    a = mod(a, m);
    if (a == 0) throw std::invalid_argument("Inverse of 0");
    long long res = 1;
    long long base = a;
    int exp = m - 2;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % m;
        base = (base * base) % m;
        exp >>= 1;
    }
    return static_cast<int>(res);
}

std::vector<int> mul_by_x_minus_a(const std::vector<int>& poly, int a, int field) {
    unsigned d = poly.size();
    std::vector<int> out(d + 1);
    for (unsigned j = 0; j <= d; ++j) out[j] = 0;

    for (unsigned j = 0; j <= d; ++j) {
        long long val = 0;
        if (j > 0) val += poly[j - 1];
        if (j < d) {
            val -= 1LL * a * poly[j];
        }
        val %= field;
        if (val < 0) val += field;
        out[j] = static_cast<int>(val);
    }
    return out;
}

// Reconstruct the secret from k shares using Lagrange interpolation
data_t reconstructFromShares(const std::vector<data_t>& shares, const std::vector<int>& x_values, unsigned k, unsigned field, unsigned kn = 1) {
    if (shares.size() < k || x_values.size() < k) {
        throw std::invalid_argument("Not enough shares or x-values to reconstruct the secret");
    }

    unsigned shareSize = shares[0].size();
    data_t result;
    result.reserve(shareSize * kn);

    // Use the provided x-values
    std::vector<int> xs(x_values.begin(), x_values.begin() + k);

    int p = static_cast<int>(field);

    for (unsigned block = 0; block < shareSize; ++block) {
        std::vector<int> coeffs(k, 0);

        for (unsigned m = 0; m < k; ++m) {
            int x_m = xs[m];
            int y_m = static_cast<int>(shares[m][block]) % p;
            if (y_m < 0) y_m += p;

            std::vector<int> poly; poly.push_back(1);

            int denom = 1;
            for (unsigned l = 0; l < k; ++l) {
                if (l == m) continue;
                int x_l = xs[l];

                poly = mul_by_x_minus_a(poly, x_l, p);

                denom = static_cast<int>((1LL * denom * mod(x_m - x_l, p)) % p);
            }

            if (denom == 0) {
                throw std::runtime_error("Singular denominator in Lagrange basis (duplicate x?)");
            }

            int invDenom = modInverseInt(denom, p); // Replaces division
            long long scale = (1LL * y_m * invDenom) % p;
            for (unsigned j = 0; j < poly.size(); ++j) {
                long long add = (scale * poly[j]) % p;
                coeffs[j] = mod(coeffs[j] + static_cast<int>(add), p);
            }
        }

        for (unsigned j = 0; j < kn; ++j) {
            int v = coeffs[j] % p;
            if (v < 0) v += p;
            result.push_back(static_cast<std::uint8_t>(v));
        }
    }

    return result;
}

int main() {
    int width, height, channels;
    std::uint8_t* img = stbi_load("input.png", &width, &height, &channels, 0);
    if (!img) {
        std::cerr << "Failed to load image\n";
        return 1;
    }
    std::size_t size = static_cast<std::size_t>(width) * height * channels;
    std::vector<std::uint8_t> imgVector(img, img + size);


    if (!img) {
        std::cerr << "Failed to load image\n";
        return 1;
    }

    unsigned k = 4;
    unsigned n = 8;
    unsigned kn = 1;
    unsigned field = 251;

    unsigned shareWidth = width;
    unsigned shareHeight = height / kn;

    std::vector<data_t> shares = getShares(imgVector, k, n, field, kn);

    for (unsigned i = 0; i < shares.size(); ++i) {
        std::string filename = "share_" + std::to_string(i + 1) + ".png";
        if (!stbi_write_png(filename.c_str(), shareWidth, shareHeight, channels, shares[i].data(), shareWidth * channels)) {
            std::cerr << "Failed to write share image: " << filename << "\n";
            return 1;
        }
    }

    std::vector<int> x_values1 = {1, 2, 3, 4};
    std::vector<data_t> sharesToUse1 = {shares[0], shares[1], shares[2], shares[3]};
    data_t reconstructed = reconstructFromShares(sharesToUse1, x_values1, k, field, kn);

    std::vector<int> x_values2 = {2, 4, 5, 6};
    std::vector<data_t> sharesToUse2 = { shares[1], shares[3], shares[4], shares[5] };
    data_t reconstructed2 = reconstructFromShares(sharesToUse2, x_values2, k, field, kn);

    if (!stbi_write_png("reconstructed.png", width, height, channels, reconstructed.data(), width * channels)) {
        std::cerr << "Failed to write reconstructed image\n";
        return 1;
    }

    if (!stbi_write_png("reconstructed_2.png", width, height, channels, reconstructed2.data(), width * channels)) {
        std::cerr << "Failed to write reconstructed image\n";
        return 1;
    }

    double psnr = computePSNR(imgVector, reconstructed);
    std::cout << "PSNR = " << psnr << " dB\n";

    psnr = computePSNR(imgVector, reconstructed2);
    std::cout << "PSNR = " << psnr << " dB\n";

    if (kn == 1) {
        std::vector<int> x_values3 = {2, 5, 6};
        std::vector<data_t> sharesToUse3 = { shares[1], shares[4], shares[5] };
        data_t reconstructed3 = reconstructFromShares(sharesToUse3, x_values2, k - 1, field, kn);

        if (!stbi_write_png("reconstructed_3.png", width, height, channels, reconstructed3.data(), width * channels)) {
            std::cerr << "Failed to write reconstructed image\n";
            return 1;
        }

        psnr = computePSNR(imgVector, reconstructed3);
        std::cout << "PSNR = " << psnr << " dB\n";
    }

    data_t imageData(img, img + width * height * channels);
    stbi_image_free(img);
}