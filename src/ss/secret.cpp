#include "ss/secret.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <limits>
#include <random>

namespace ss {

// Field helper implementations

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

// PrimeField

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

// BinaryField

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
typename BinaryField<Storage>::storage_type BinaryField<Storage>::sub(storage_type a, storage_type b) const noexcept {
    return add(a, b);
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

// getShares helper & wrappers

template<typename Storage, typename CoeffGenerator>
static std::vector<std::vector<Storage>> getSharesHelper(
    const std::vector<Storage>& data,
    unsigned k,
    unsigned n,
    const Field<Storage>& F,
    unsigned kn,
    CoeffGenerator coeffGenerator,
    const std::vector<Storage>& evalPoints
) {
    if (k == 0 || n == 0) throw std::invalid_argument("k and n must be > 0");
    if (kn == 0 || kn > k) throw std::invalid_argument("kn must satisfy 1 <= kn <= k");

    if (F.is_binary_field()) {
        const auto& BF = static_cast<const BinaryField<Storage>&>(F);
        if (n > ((1ULL << BF.m) - 1ULL))
            throw std::invalid_argument("n too large for field size");
    }

    if (!evalPoints.empty()) {
        if (evalPoints.size() != n)
            throw std::invalid_argument("evalPoints size must equal n");

        std::unordered_set<std::uint64_t> seen;
        std::uint64_t order = F.getOrder();
        for (auto x : evalPoints) {
            std::uint64_t xv = static_cast<std::uint64_t>(x);
            if (xv == 0) throw std::invalid_argument("evaluation points must be non-zero");
            if (xv >= order) throw std::invalid_argument("evaluation point out of field range");
            if (!seen.insert(xv).second)
                throw std::invalid_argument("evaluation points must be distinct");
        }
    }

    const unsigned shareSize = static_cast<unsigned>((data.size() + kn - 1) / kn);
    std::vector<std::vector<Storage>> shares(n);
    for (auto &v : shares) v.reserve(shareSize);

    std::random_device dev;
    std::mt19937_64 rng(dev());

    std::uint64_t max_val;
    if (F.is_binary_field()) {
        const auto& BF = static_cast<const BinaryField<Storage>&>(F);
        max_val = (BF.m >= 64) ? std::numeric_limits<std::uint64_t>::max() : ((1ULL << BF.m) - 1ULL);
    } else {
        const auto& PF = static_cast<const PrimeField<Storage>&>(F);
        max_val = PF.p - 1;
    }
    std::uniform_int_distribution<std::uint64_t> dist(0, max_val);

    unsigned dataCursor = 0;
    for (unsigned block = 0; block < shareSize; ++block) {
        auto polyCoeffs = coeffGenerator(k, kn, rng, dist);

        // Pack secret data into coefficients
        for (unsigned j = 0; j < kn; ++j) {
            if (dataCursor < data.size()) {
                Storage value = data[dataCursor++];

                // Reduce mod field prime
                if (!F.is_binary_field()) {
                    const auto& PF = static_cast<const PrimeField<Storage>&>(F);
                    std::uint64_t v = static_cast<std::uint64_t>(value);
                    if (v >= PF.p) v %= PF.p;
                    value = static_cast<Storage>(v);
                }
                polyCoeffs[j] = value;
            } else {
                polyCoeffs[j] = static_cast<Storage>(0);
            }
        }

        // Evaluate polynomial
        for (unsigned j = 0; j < n; ++j) {
            Storage x = evalPoints.empty() ? static_cast<Storage>(j + 1) : evalPoints[j];
            Storage y = F.evaluatePolynomial(polyCoeffs, x);
            shares[j].push_back(y);
        }
    }
    return shares;
}

template<typename Storage>
std::vector<std::vector<Storage>> getShares(
    const std::vector<Storage>& data,
    unsigned k,
    unsigned n,
    const Field<Storage>& F,
    unsigned kn,
    const std::vector<Storage>& evalPoints
) {
    auto coeffGenerator = [](unsigned k, unsigned kn, std::mt19937_64& rng, std::uniform_int_distribution<std::uint64_t>& dist) {
        std::vector<Storage> polyCoeffs(k);
        // Fill all coefficients with random values
        for (unsigned j = 0; j < k; ++j) {
            polyCoeffs[j] = static_cast<Storage>(dist(rng));
        }
        return polyCoeffs;
    };

    return getSharesHelper(data, k, n, F, kn, coeffGenerator, evalPoints);
}

template<typename Storage>
std::vector<std::vector<Storage>> getSharesSymmetric(
    const std::vector<Storage>& data,
    unsigned k,
    unsigned n,
    const Field<Storage>& F,
    Storage a,
    unsigned kn,
    const std::vector<Storage>& evalPoints
) {
    if (n % 2 == 0) throw std::invalid_argument("n must be odd for symmetric shares");
    if (k % 2 == 0) throw std::invalid_argument("k must be odd for symmetric shares");

    auto binomial_u64 = [](unsigned n, unsigned r) -> std::uint64_t {
        if (r > n) return 0;
        unsigned kmin = (r > n - r) ? n - r : r;
        std::uint64_t res = 1;
        for (unsigned i = 1; i <= kmin; ++i) {
            res = res * (n - kmin + i) / i;
        }
        return res;
    };

    auto coeffGenerator = [a, &binomial_u64](unsigned k, unsigned /*kn*/, std::mt19937_64& rng, std::uniform_int_distribution<std::uint64_t>& dist) {
        std::vector<Storage> polyCoeffs(k, static_cast<Storage>(0));

        unsigned num_coeffs = (k + 1) / 2;

        std::vector<Storage> c(num_coeffs, static_cast<Storage>(0));
        for (unsigned j = 0; j < num_coeffs; ++j) {
            c[j] = static_cast<Storage>(dist(rng));
        }

        for (unsigned j = 0; j < num_coeffs; ++j) {
            unsigned power = 2 * j;
            Storage cj = c[j];

            for (unsigned r = 0; r <= power && r < k; ++r) {
                std::uint64_t bin = binomial_u64(power, r);
                Storage term = cj * static_cast<Storage>(bin);

                unsigned exp = power - r;
                Storage a_pow = static_cast<Storage>(1);
                for (unsigned e = 0; e < exp; ++e) {
                    a_pow = a_pow * a;
                }
                if ((exp & 1u) == 1u) {
                    a_pow = static_cast<Storage>(0) - a_pow;
                }

                term = term * a_pow;

                polyCoeffs[r] = polyCoeffs[r] + term;
            }
        }

        return polyCoeffs;
    };

    return getSharesHelper(data, k, n, F, kn, coeffGenerator, evalPoints);
}

template<typename Storage>
std::vector<Storage> lagrangeBasisCoeffs(
    const std::vector<Storage>& xs,
    unsigned index,
    const Field<Storage>& F
) {
    if (xs.empty()) throw std::invalid_argument("xs must be non-empty");
    if (index >= xs.size()) throw std::out_of_range("index out of range");

    std::vector<Storage> poly(1, static_cast<Storage>(1));
    for (unsigned l = 0; l < xs.size(); ++l) {
        if (l == index) continue;
        poly = F.mul_by_x_minus_a(poly, xs[l]);
    }

    Storage denom = static_cast<Storage>(1);
    for (unsigned l = 0; l < xs.size(); ++l) {
        if (l == index) continue;
        Storage diff = F.sub(xs[index], xs[l]);
        denom = F.mul(denom, diff);
    }

    if (denom == static_cast<Storage>(0)) {
        throw std::runtime_error("Singular denominator in Lagrange basis (duplicate x?)");
    }

    Storage invDenom = F.inv(denom);

    for (size_t j = 0; j < poly.size(); ++j) {
        poly[j] = F.mul(poly[j], invDenom);
    }

    return poly;
}

template<typename Storage>
Storage evaluateLagrangeBasisAt(
    const std::vector<Storage>& xs,
    unsigned index,
    Storage x,
    const Field<Storage>& F
) {
    return F.evaluatePolynomial(lagrangeBasisCoeffs<Storage>(xs, index, F), x);
}

template<typename Storage>
std::vector<Storage> reconstructFromShares(
    const std::vector<std::vector<Storage>>& shares,
    const std::vector<Storage>& evalPoints,
    unsigned k,
    const Field<Storage>& F,
    unsigned kn,
    std::size_t expected_size
) {
    if (shares.size() < k) throw std::invalid_argument("Not enough shares");
    if (k == 0) throw std::invalid_argument("k must be > 0");

    if (evalPoints.size() < k) throw std::invalid_argument("Not enough evalPoints");

    const unsigned shareSize = static_cast<unsigned>(shares[0].size());
    for (unsigned i = 0; i < k; ++i) {
        if (shares[i].size() != shareSize) {
            throw std::invalid_argument("Share lengths differ");
        }
    }

    std::vector<Storage> result;
    result.reserve(static_cast<size_t>(shareSize) * kn);

    // Use first k evaluation points
    std::vector<Storage> xs(k);
    for (unsigned i = 0; i < k; ++i) {
        xs[i] = evalPoints[i];
    }

    for (unsigned block = 0; block < shareSize; ++block) {
        std::vector<Storage> coeffs(k, static_cast<Storage>(0));

        for (unsigned index = 0; index < k; ++index) {
            Storage y_m = shares[index][block];

            auto basis = lagrangeBasisCoeffs<Storage>(xs, index, F);

            // Add y_m * basis to coeffs
            for (size_t j = 0; j < basis.size(); ++j) {
                Storage addv = F.mul(y_m, basis[j]);
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

template<typename Storage>
std::pair<std::vector<std::vector<Storage>>, std::vector<Storage>> selectSharesAndEvalPoints(
    const std::vector<unsigned>& indices,
    const std::vector<std::vector<Storage>>& shares,
    const std::vector<Storage>& evalPoints
) {
    if (indices.size() > shares.size()) {
        throw std::invalid_argument("Too many indices");
    }

    std::vector<Storage> actual_evalPoints;
    if (evalPoints.empty()) {
        actual_evalPoints.resize(shares.size());
        std::iota(actual_evalPoints.begin(), actual_evalPoints.end(), static_cast<Storage>(1));
    } else {
        actual_evalPoints = evalPoints;
    }

    std::vector<std::vector<Storage>> selected_shares;
    std::vector<Storage> selected_evalPoints;
    selected_shares.reserve(indices.size());
    selected_evalPoints.reserve(indices.size());
    for (unsigned idx : indices) {
        if (idx >= shares.size()) {
            throw std::out_of_range("Index out of range in selectSharesAndEvalPoints");
        }
        selected_shares.push_back(shares[idx]);
        if (idx >= actual_evalPoints.size()) {
            throw std::out_of_range("Index out of range in evalPoints in selectSharesAndEvalPoints");
        }
        selected_evalPoints.push_back(actual_evalPoints[idx]);
    }
    return {selected_shares, selected_evalPoints};
}

// Explicit template instantiations

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
template std::vector<std::vector<std::uint8_t>> getShares<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, unsigned, const Field<std::uint8_t>&, unsigned, const std::vector<std::uint8_t>&);
template std::vector<std::vector<std::uint16_t>> getShares<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, unsigned, const Field<std::uint16_t>&, unsigned, const std::vector<std::uint16_t>&);
template std::vector<std::vector<std::uint32_t>> getShares<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, unsigned, const Field<std::uint32_t>&, unsigned, const std::vector<std::uint32_t>&);

template std::vector<std::vector<std::uint8_t>> getSharesSymmetric<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, unsigned, const Field<std::uint8_t>&, std::uint8_t a, unsigned, const std::vector<std::uint8_t>&);
template std::vector<std::vector<std::uint16_t>> getSharesSymmetric<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, unsigned, const Field<std::uint16_t>&, std::uint16_t a, unsigned, const std::vector<std::uint16_t>&);
template std::vector<std::vector<std::uint32_t>> getSharesSymmetric<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, unsigned, const Field<std::uint32_t>&, std::uint32_t a, unsigned, const std::vector<std::uint32_t>&);

template std::vector<std::uint8_t> reconstructFromShares<std::uint8_t>(const std::vector<std::vector<std::uint8_t>>&, const std::vector<std::uint8_t>&, unsigned, const Field<std::uint8_t>&, unsigned, std::size_t);
template std::vector<std::uint16_t> reconstructFromShares<std::uint16_t>(const std::vector<std::vector<std::uint16_t>>&, const std::vector<std::uint16_t>&, unsigned, const Field<std::uint16_t>&, unsigned, std::size_t);
template std::vector<std::uint32_t> reconstructFromShares<std::uint32_t>(const std::vector<std::vector<std::uint32_t>>&, const std::vector<std::uint32_t>&, unsigned, const Field<std::uint32_t>&, unsigned, std::size_t);

template std::pair<std::vector<std::vector<std::uint8_t>>, std::vector<std::uint8_t>> selectSharesAndEvalPoints(const std::vector<unsigned>& indices, const std::vector<std::vector<std::uint8_t>>& shares, const std::vector<std::uint8_t>& evalPoints);
template std::pair<std::vector<std::vector<std::uint16_t>>, std::vector<std::uint16_t>> selectSharesAndEvalPoints(const std::vector<unsigned>& indices, const std::vector<std::vector<std::uint16_t>>& shares, const std::vector<std::uint16_t>& evalPoints);
template std::pair<std::vector<std::vector<std::uint32_t>>, std::vector<std::uint32_t>> selectSharesAndEvalPoints(const std::vector<unsigned>& indices, const std::vector<std::vector<std::uint32_t>>& shares, const std::vector<std::uint32_t>& evalPoints);

template std::vector<std::uint8_t> lagrangeBasisCoeffs<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, const Field<std::uint8_t>&);
template std::vector<std::uint16_t> lagrangeBasisCoeffs<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, const Field<std::uint16_t>&);
template std::vector<std::uint32_t> lagrangeBasisCoeffs<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, const Field<std::uint32_t>&);

template std::uint8_t evaluateLagrangeBasisAt<std::uint8_t>(const std::vector<std::uint8_t>&, unsigned, std::uint8_t, const Field<std::uint8_t>&);
template std::uint16_t evaluateLagrangeBasisAt<std::uint16_t>(const std::vector<std::uint16_t>&, unsigned, std::uint16_t, const Field<std::uint16_t>&);
template std::uint32_t evaluateLagrangeBasisAt<std::uint32_t>(const std::vector<std::uint32_t>&, unsigned, std::uint32_t, const Field<std::uint32_t>&);

} // namespace ss