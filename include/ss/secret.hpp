#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <numeric>

namespace ss {

static constexpr std::uint64_t irr[32] = {
    0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83, 0x11b, 0x203, 0x409, 0x805, 0x1009, 0x201b, 0x4021, 0x8003, 0x1002b, 0x20003, 0x40009, 0x80027, 0x100009, 0x200005, 0x400003, 0x800021, 0x100001b, 0x2000009, 0x400001b, 0x8000027, 0x10000003, 0x20000005, 0x40000003, 0x80000009, 0x10000008d
};

// Galois Field interface
template<typename Storage>
struct Field {
    static_assert(std::is_unsigned<Storage>::value, "Storage must be unsigned integer type");
    using storage_type = Storage;

    virtual ~Field() = default;
    
    virtual storage_type add(storage_type a, storage_type b) const noexcept = 0;
    virtual storage_type sub(storage_type a, storage_type b) const noexcept = 0;
    virtual storage_type mul(storage_type a, storage_type b) const noexcept = 0;
    virtual storage_type pow(storage_type a, std::uint64_t e) const = 0;
    virtual storage_type inv(storage_type a) const = 0;

    virtual std::uint64_t getOrder() const = 0;
    
    storage_type evaluatePolynomial(const std::vector<storage_type>& coefficients, storage_type x) const;
    std::vector<storage_type> mul_by_x_minus_a(const std::vector<storage_type>& poly, storage_type a) const;
    
    virtual bool is_binary_field() const noexcept = 0;
};

// Prime field GF(p)
template<typename Storage>
struct PrimeField : public Field<Storage> {
    using typename Field<Storage>::storage_type;

    std::uint64_t p;

    explicit PrimeField(std::uint64_t prime);
    
    storage_type add(storage_type a, storage_type b) const noexcept override;
    storage_type sub(storage_type a, storage_type b) const noexcept override;
    storage_type mul(storage_type a, storage_type b) const noexcept override;
    storage_type pow(storage_type a, std::uint64_t e) const override;
    storage_type inv(storage_type a) const override;

    std::uint64_t getOrder() const override { return p; }
    
    bool is_binary_field() const noexcept override { return false; }
};

// Binary extension field GF(2^m)
template<typename Storage>
struct BinaryField : public Field<Storage> {
    using typename Field<Storage>::storage_type;
    unsigned m;

    BinaryField(unsigned m_deg);
    
    storage_type add(storage_type a, storage_type b) const noexcept override;
    storage_type sub(storage_type a, storage_type b) const noexcept override;
    storage_type mul(storage_type a, storage_type b) const noexcept override;
    storage_type pow(storage_type a, std::uint64_t e) const override;
    storage_type inv(storage_type a) const override;

    std::uint64_t getOrder() const override { return static_cast<std::uint64_t>(1) << m; }
    
    bool is_binary_field() const noexcept override { return true; }
};

template<typename Storage>
std::vector<Storage> lagrangeBasisCoeffs(const std::vector<Storage>& xs, unsigned index, const Field<Storage>& F);

template<typename Storage>
Storage evaluateLagrangeBasisAt(const std::vector<Storage>& xs, unsigned index, Storage x, const Field<Storage>& F);

template<typename Storage>
std::vector<std::vector<Storage>> getShares(const std::vector<Storage>& data, unsigned k, unsigned n, const Field<Storage>& F, unsigned kn = 1, const std::vector<Storage>& evalPoints = {});

template<typename Storage>
std::vector<std::vector<Storage>> getSharesSymmetric(const std::vector<Storage>& data, unsigned k, unsigned n, const Field<Storage>& F, Storage a, unsigned kn = 1, const std::vector<Storage>& evalPoints = {});

template<typename Storage>
std::vector<Storage> reconstructFromShares(const std::vector<std::vector<Storage>>& shares, const std::vector<Storage>& evalPoints, unsigned k, const Field<Storage>& F, unsigned kn = 1, std::size_t expected_size = static_cast<std::size_t>(-1));

template<typename Storage>
std::pair<std::vector<std::vector<Storage>>, std::vector<Storage>> selectSharesAndEvalPoints(const std::vector<unsigned>& indices, const std::vector<std::vector<Storage>>& shares, const std::vector<Storage>& evalPoints = {});

} // namespace ss