#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace ss {

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
    storage_type sub(storage_type a, storage_type b) const noexcept override { return add(a, b); }
    storage_type mul(storage_type a, storage_type b) const noexcept override;
    storage_type pow(storage_type a, std::uint64_t e) const override;
    storage_type inv(storage_type a) const override;

    std::uint64_t getOrder() const override { return static_cast<std::uint64_t>(1) << m; }
    
    bool is_binary_field() const noexcept override { return true; }
};

template<typename Storage>
std::vector<std::vector<Storage>> getShares(const std::vector<Storage>& data, unsigned k, unsigned n, const Field<Storage>& F, unsigned kn = 1, bool clampToOrder = true);

template<typename Storage>
std::vector<Storage> reconstructFromShares(const std::vector<std::vector<Storage>>& shares, const std::vector<Storage>& x_values, unsigned k, const Field<Storage>& F, unsigned kn = 1, std::size_t expected_size = static_cast<std::size_t>(-1));

} // namespace ss