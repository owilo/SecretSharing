#include "secret.hpp"

#include <vector>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <random>
#include <type_traits>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ss {

static constexpr std::uint64_t irr[32] = {
    0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83, 0x11b, 0x203, 0x409, 0x805, 0x1009, 0x201b, 0x4021, 0x8003, 0x1002b, 0x20003, 0x40009, 0x80027, 0x100009, 0x200005, 0x400003, 0x800021, 0x100001b, 0x2000009, 0x400001b, 0x8000027, 0x10000003, 0x20000005, 0x40000003, 0x80000009, 0x10000008d
};

// Field definitions

template<typename Storage>
typename Field<Storage>::storage_type Field<Storage>::evaluatePolynomial(const std::vector<storage_type>& coefficients, storage_type x) const {
    storage_type acc = static_cast<storage_type>(0);
    for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
        acc = add(mul(acc, x), coefficients[i]);
    }
    return acc;
}

template<typename Storage>
std::vector<typename Field<Storage>::storage_type> Field<Storage>::mul_by_x_minus_a(const std::vector<storage_type>& poly, storage_type a) const {
    const std::size_t d = poly.size();
    std::vector<storage_type> out(d + 1, static_cast<storage_type>(0));
    for (std::size_t j = 0; j <= d; ++j) {
        storage_type val = static_cast<storage_type>(0);
        if (j > 0) val = add(val, poly[j - 1]);
        if (j < d) val = sub(val, mul(a, poly[j]));
        out[j] = val;
    }
    return out;
}

// PrimeField definitions

template<typename Storage>
PrimeField<Storage>::PrimeField(std::uint64_t prime) : p(prime) {
    if (p < 2) throw std::invalid_argument("Prime must be >= 2");
}

template<typename Storage>
typename PrimeField<Storage>::storage_type PrimeField<Storage>::add(storage_type a, storage_type b) const noexcept {
    std::uint64_t s = static_cast<std::uint64_t>(a) + static_cast<std::uint64_t>(b);
    return static_cast<storage_type>(s % p);
}

template<typename Storage>
typename PrimeField<Storage>::storage_type PrimeField<Storage>::sub(storage_type a, storage_type b) const noexcept {
    std::int64_t v = static_cast<std::int64_t>(a) - static_cast<std::int64_t>(b);
    std::int64_t r = v % static_cast<std::int64_t>(p);
    if (r < 0) r += static_cast<std::int64_t>(p);
    return static_cast<storage_type>(r);
}

template<typename Storage>
typename PrimeField<Storage>::storage_type PrimeField<Storage>::mul(storage_type a, storage_type b) const noexcept {
    std::uint64_t r = static_cast<std::uint64_t>(a) * static_cast<std::uint64_t>(b);
    return static_cast<storage_type>(r % p);
}

template<typename Storage>
typename PrimeField<Storage>::storage_type PrimeField<Storage>::pow(storage_type a, std::uint64_t e) const {
    if (p == 2) return a & 1;
    std::uint64_t base = static_cast<std::uint64_t>(a) % p;
    std::uint64_t res = 1 % p;
    while (e) {
        if (e & 1) res = (res * base) % p;
        base = (base * base) % p;
        e >>= 1;
    }
    return static_cast<storage_type>(res);
}

template<typename Storage>
typename PrimeField<Storage>::storage_type PrimeField<Storage>::inv(storage_type a) const {
    if (a == 0) throw std::invalid_argument("PrimeField::inv of zero");
    return pow(a, p - 2);
}

// BinaryField definitions

template<typename Storage>
BinaryField<Storage>::BinaryField(unsigned m_deg)
    : m(m_deg) {
    if (m == 0 || m > (8 * sizeof(Storage))) throw std::invalid_argument("Invalid extension degree m");
}

template<typename Storage>
typename BinaryField<Storage>::storage_type BinaryField<Storage>::add(storage_type a, storage_type b) const noexcept {
    return static_cast<storage_type>(a ^ b);
}

template<typename Storage>
typename BinaryField<Storage>::storage_type BinaryField<Storage>::mul(storage_type a, storage_type b) const noexcept {
    using U = std::uint64_t;
    U aa = static_cast<U>(a);
    U bb = static_cast<U>(b);
    U res = 0;
    while (bb) {
        if (bb & 1u) res ^= aa;
        bb >>= 1;
        aa <<= 1;
    }
    // Reduce modulo the polynomial
    for (int bit = static_cast<int>(2 * m - 2); bit >= static_cast<int>(m); --bit) {
        if (res & (static_cast<U>(1) << bit)) {
            res ^= (static_cast<U>(irr[m - 1]) << (bit - static_cast<int>(m)));
        }
    }
    U mask = ((static_cast<U>(1) << m) - 1u);
    return static_cast<storage_type>(res & mask);
}

template<typename Storage>
typename BinaryField<Storage>::storage_type BinaryField<Storage>::pow(storage_type a, std::uint64_t e) const {
    storage_type r = static_cast<storage_type>(1);
    while (e) {
        if (e & 1) r = mul(r, a);
        a = mul(a, a);
        e >>= 1;
    }
    return r;
}

template<typename Storage>
typename BinaryField<Storage>::storage_type BinaryField<Storage>::inv(storage_type a) const {
    if (a == 0) throw std::invalid_argument("BinaryField::inv of zero");
    if (m >= 64) throw std::runtime_error("m too large to compute inverse via exponentiation safely");
    std::uint64_t e = (static_cast<std::uint64_t>(1) << m) - 2u;
    return pow(a, e);
}

// Secret sharing

template<typename Storage, typename CoeffGenerator>
std::vector<std::vector<Storage>> getSharesHelper(
    const std::vector<Storage>& data, 
    unsigned k, 
    unsigned n, 
    const Field<Storage>& F, 
    unsigned kn, 
    bool clampToOrder,
    CoeffGenerator coeffGenerator
) {
    if (k == 0 || n == 0) throw std::invalid_argument("k and n must be > 0");
    if (kn == 0 || kn > k) throw std::invalid_argument("kn must satisfy 1 <= kn <= k");
    
    if (F.is_binary_field()) {
        const auto& BF = static_cast<const BinaryField<Storage>&>(F);
        if (n > ((1ULL << BF.m) - 1ULL)) throw std::invalid_argument("n too large for field size");
    }

    const unsigned shareSize = static_cast<unsigned>((data.size() + kn - 1) / kn);
    std::vector<std::vector<Storage>> shares(n);
    for (unsigned i = 0; i < n; ++i) shares[i].reserve(shareSize);

    std::random_device dev;
    std::mt19937_64 rng(dev());
    
    std::uint64_t max_val;
    if (F.is_binary_field()) {
        const auto& BF = static_cast<const BinaryField<Storage>&>(F);
        if (BF.m >= 64) {
            max_val = std::numeric_limits<std::uint64_t>::max();
        } else {
            max_val = (1ULL << BF.m) - 1ULL;
        }
    } else {
        const auto& PF = static_cast<const PrimeField<Storage>&>(F);
        max_val = PF.p - 1;
    }

    std::uniform_int_distribution<std::uint64_t> dist(0, max_val);

    unsigned dataCursor = 0;
    for (unsigned block = 0; block < shareSize; ++block) {
        // Get polynomial coefficients from the provided generator
        std::vector<Storage> polyCoeffs = coeffGenerator(k, kn, rng, dist);
        
        // Pack secret data into the first kn coefficients
        for (unsigned j = 0; j < kn; ++j) {
            if (dataCursor < data.size()) {
                Storage value = data[dataCursor++];
                if (clampToOrder) {
                    auto order = F.getOrder();
                    if (static_cast<std::uint64_t>(value) >= order) {
                        value = static_cast<Storage>(order - 1);
                    }
                } else if (!F.is_binary_field()) {
                    const auto& PF = static_cast<const PrimeField<Storage>&>(F);
                    auto v = static_cast<std::uint64_t>(value);
                    if (v >= PF.p) v = v % PF.p;
                    value = static_cast<Storage>(v);
                }
                polyCoeffs[j] = value;
            } else {
                polyCoeffs[j] = static_cast<Storage>(0);
            }
        }

        // Evaluate polynomial
        for (unsigned j = 0; j < n; ++j) {
            Storage x = static_cast<Storage>(j + 1);
            Storage y = F.evaluatePolynomial(polyCoeffs, x);
            shares[j].push_back(y);
        }
    }
    return shares;
}

template<typename Storage>
std::vector<std::vector<Storage>> getShares(const std::vector<Storage>& data, unsigned k, unsigned n, const Field<Storage>& F, unsigned kn, bool clampToOrder) {
    auto coeffGenerator = [](unsigned k, unsigned kn, std::mt19937_64& rng, std::uniform_int_distribution<std::uint64_t>& dist) {
        std::vector<Storage> polyCoeffs(k);
        // Fill all coefficients with random values
        for (unsigned j = 0; j < k; ++j) {
            polyCoeffs[j] = static_cast<Storage>(dist(rng));
        }
        return polyCoeffs;
    };
    
    return getSharesHelper(data, k, n, F, kn, clampToOrder, coeffGenerator);
}

template<typename Storage>
std::vector<std::vector<Storage>> getSharesSymmetric(const std::vector<Storage>& data, unsigned k, unsigned n, const Field<Storage>& F, unsigned kn, bool clampToOrder) {
    if (n % 2 == 0) throw std::invalid_argument("n must be odd for symmetric shares");
    if (k % 2 == 0) throw std::invalid_argument("k must be odd for symmetric shares");
    
    auto coeffGenerator = [](unsigned k, unsigned kn, std::mt19937_64& rng, std::uniform_int_distribution<std::uint64_t>& dist) {
        std::vector<Storage> polyCoeffs(k, static_cast<Storage>(0));
        
        unsigned num_coeffs = (k + 1) / 2;
        
        // Generate random coefficients for even powers only
        for (unsigned j = 0; j < num_coeffs; ++j) {
            unsigned power = 2 * j;
            if (power < k) {
                polyCoeffs[power] = static_cast<Storage>(dist(rng));
            }
        }
        
        return polyCoeffs;
    };
    
    return getSharesHelper(data, k, n, F, kn, clampToOrder, coeffGenerator);
}

template<typename Storage>
std::vector<Storage> reconstructFromShares(const std::vector<std::vector<Storage>>& shares, const std::vector<Storage>& x_values, unsigned k, const Field<Storage>& F, unsigned kn, std::size_t expected_size) {
    if (shares.size() < k) throw std::invalid_argument("Not enough shares");
    if (x_values.size() < k) throw std::invalid_argument("Not enough x_values");
    
    const unsigned shareSize = static_cast<unsigned>(shares[0].size());
    for (unsigned i = 0; i < k; ++i) {
        if (shares[i].size() != shareSize) {
            throw std::invalid_argument("Share lengths differ");
        }
    }

    std::vector<Storage> result;
    result.reserve(static_cast<size_t>(shareSize) * kn);

    std::vector<Storage> xs(k);
    for (unsigned i = 0; i < k; ++i) {
        xs[i] = x_values[i];
    }

    for (unsigned block = 0; block < shareSize; ++block) {
        std::vector<Storage> coeffs(k, static_cast<Storage>(0));

        for (unsigned m_idx = 0; m_idx < k; ++m_idx) {
            Storage x_m = xs[m_idx];
            Storage y_m = shares[m_idx][block];
            
            std::vector<Storage> poly(1, static_cast<Storage>(1));
            Storage denom = static_cast<Storage>(1);
            
            // Build Lagrange basis polynomial and denominator
            for (unsigned l = 0; l < k; ++l) {
                if (l == m_idx) continue;
                poly = F.mul_by_x_minus_a(poly, xs[l]);
                Storage diff = F.sub(x_m, xs[l]);
                denom = F.mul(denom, diff);
            }
            
            if (denom == 0) throw std::runtime_error("Singular denominator in Lagrange basis (duplicate x?)");
            
            Storage invDenom = F.inv(denom);
            Storage scale = F.mul(y_m, invDenom);
            
            // Accumulate the scaled basis polynomial
            for (size_t j = 0; j < poly.size(); ++j) {
                Storage addv = F.mul(scale, poly[j]);
                coeffs[j] = F.add(coeffs[j], addv);
            }
        }

        // Extract the secret coefficients
        for (unsigned j = 0; j < kn; ++j) {
            result.push_back(coeffs[j]);
        }
    }

    const std::size_t AUTO = static_cast<std::size_t>(-1);
    if (expected_size == AUTO) {
        while (!result.empty() && result.back() == static_cast<Storage>(0)) {
            result.pop_back();
        }
    } else {
        if (result.size() > expected_size) {
            result.resize(expected_size);
        } else if (result.size() < expected_size) {
            result.resize(expected_size, static_cast<Storage>(0));
        }
    }

    return result;
}

// Template instantiations

// Field base class instantiations
template struct Field<std::uint8_t>;
template struct Field<std::uint16_t>;
template struct Field<std::uint32_t>;

// PrimeField class instantiations
template class PrimeField<std::uint8_t>;
template class PrimeField<std::uint16_t>;
template class PrimeField<std::uint32_t>;

// BinaryField class instantiations
template class BinaryField<std::uint8_t>;
template class BinaryField<std::uint16_t>;
template class BinaryField<std::uint32_t>;

// Secret sharing function instantiations
template std::vector<std::vector<std::uint8_t>> getShares<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, unsigned, const Field<std::uint8_t>&, unsigned, bool);
template std::vector<std::vector<std::uint16_t>> getShares<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, unsigned, const Field<std::uint16_t>&, unsigned, bool);
template std::vector<std::vector<std::uint32_t>> getShares<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, unsigned, const Field<std::uint32_t>&, unsigned, bool);

template std::vector<std::vector<std::uint8_t>> getSharesSymmetric<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, unsigned, const Field<std::uint8_t>&, unsigned, bool);
template std::vector<std::vector<std::uint16_t>> getSharesSymmetric<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, unsigned, const Field<std::uint16_t>&, unsigned, bool);
template std::vector<std::vector<std::uint32_t>> getSharesSymmetric<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, unsigned, const Field<std::uint32_t>&, unsigned, bool);

template std::vector<std::uint8_t> reconstructFromShares<std::uint8_t>(const std::vector<std::vector<std::uint8_t>>&, const std::vector<std::uint8_t>&, unsigned, const Field<std::uint8_t>&, unsigned, std::size_t);
template std::vector<std::uint16_t> reconstructFromShares<std::uint16_t>(const std::vector<std::vector<std::uint16_t>>&, const std::vector<std::uint16_t>&, unsigned, const Field<std::uint16_t>&, unsigned, std::size_t);
template std::vector<std::uint32_t> reconstructFromShares<std::uint32_t>(const std::vector<std::vector<std::uint32_t>>&, const std::vector<std::uint32_t>&, unsigned, const Field<std::uint32_t>&, unsigned, std::size_t);

} // namespace ss