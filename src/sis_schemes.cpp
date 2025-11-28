#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <filesystem>
#include <tuple>
#include <unordered_set>

#include "ss/secret.hpp"
#include "ss/image.hpp"

namespace fs = std::filesystem;

enum class FieldBehavior {
    CLAMP,
    STRETCH,
    KEEP_PLAIN,
    SPLIT
};

static std::vector<std::size_t> indicesWhere(const std::vector<std::uint8_t>& vec, std::function<bool(std::uint8_t)> pred) {
    std::vector<std::size_t> idx;
    idx.reserve(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (pred(vec[i])) idx.push_back(i);
    }
    return idx;
}

static std::vector<unsigned> makeIndicesA(unsigned k) {
    std::vector<unsigned> ia(k);
    std::iota(ia.begin(), ia.end(), 0u);
    return ia;
}

static std::vector<unsigned> makeIndicesB(unsigned n, unsigned k) {
    std::vector<unsigned> idxB;
    for (unsigned i = 2; i < 2 + k && i < n; ++i) idxB.push_back(i);
    if (idxB.size() < k) {
        for (unsigned i = 0; i < n && idxB.size() < k; ++i) {
            if (std::find(idxB.begin(), idxB.end(), i) == idxB.end()) idxB.push_back(i);
        }
    }
    return idxB;
}

static void saveSharesImages(const std::string& result_prefix, const std::vector<std::vector<std::uint8_t>>& shares, int width, int height, bool divide_share_height_by_k, unsigned k) {
    fs::create_directories(result_prefix + "/shadows");
    int share_h = divide_share_height_by_k ? (height / k) : height;
    for (std::size_t i = 0; i < shares.size(); ++i) {
        std::string outp = result_prefix + "/shadows/shadow_" + std::to_string(i+1) + ".png";
        ss::saveGrayscalePNG(outp, shares[i], width, share_h);
    }
}

static void printMetrics(const std::string& path, const std::vector<std::uint8_t>& secret, const std::vector<std::uint8_t>& rec) {
    std::cout << "  " << path
              << " | PSNR : " << ss::computePSNR(secret, rec) << " dB"
              << " | NPCR : " << ss::computeNPCR(secret, rec) << "%"
              << " | UACI : " << ss::computeUACI(secret, rec) << "%"
              << " | Entropy: " << ss::computeEntropy(rec) << " bits/pixel\n";
}

static void warnShareCount(std::size_t got, std::size_t expect) {
    if (got != expect) {
        std::cerr << "Warning: getShares returned " << got << " shares (expected " << expect << ")\n";
    }
}

template <typename FieldType>
static std::vector<std::uint8_t> preprocessFieldVec(const std::vector<std::uint8_t>& secret_vec, const FieldType& field, FieldBehavior behavior) {
    const unsigned order = static_cast<unsigned>(field.getOrder());
    if (behavior == FieldBehavior::CLAMP) {
        return ss::clampPixels(secret_vec, 0, static_cast<std::uint8_t>(order - 1));
    } else if (behavior == FieldBehavior::STRETCH) {
        return ss::stretchHistogram(secret_vec, 0, 255, 0, static_cast<std::uint8_t>(order - 1));
    }
    return {};
}

template <typename FieldType>
static std::vector<std::vector<std::uint8_t>> createShares(const std::vector<std::uint8_t>& field_vec, unsigned k, unsigned n, const FieldType& field, unsigned kn) {
    auto orig_shares = ss::getShares<std::uint8_t>(field_vec, k, n, field, kn);
    warnShareCount(orig_shares.size(), n);
    return orig_shares;
}

template <typename FieldType>
static std::vector<std::vector<std::uint8_t>> createSharesKeepPlain(const FieldType& field, const std::vector<std::uint8_t>& secret_vec, unsigned n, unsigned k, unsigned kn, const std::vector<std::size_t>& in_idxs) {
    std::vector<std::uint8_t> packed_in(in_idxs.size());
    for (std::size_t i = 0; i < in_idxs.size(); ++i) packed_in[i] = secret_vec[in_idxs[i]];

    auto packed_shares = ss::getShares<std::uint8_t>(packed_in, k, n, field, kn);
    warnShareCount(packed_shares.size(), n);

    std::vector<std::vector<std::uint8_t>> shares(n, std::vector<std::uint8_t>(secret_vec.size(), 0));
    for (unsigned s = 0; s < n; ++s) {
        const auto &p = packed_shares[s];
        for (std::size_t j = 0; j < in_idxs.size(); ++j) {
            shares[s][in_idxs[j]] = p[j];
        }
    }
    return shares;
}

static std::vector<std::vector<std::uint8_t>> createSharesSplit(unsigned splitA, unsigned splitB, const std::vector<std::uint8_t>& secret_vec, unsigned n, unsigned k, unsigned kn, const std::vector<std::size_t>& idxA_in, const std::vector<std::size_t>& idxB_in) {
    ss::PrimeField<std::uint8_t> fieldA(splitA);
    ss::PrimeField<std::uint8_t> fieldB(splitB);

    std::vector<std::uint8_t> packedA(idxA_in.size()), packedB(idxB_in.size());
    for (std::size_t i = 0; i < idxA_in.size(); ++i) packedA[i] = secret_vec[idxA_in[i]];
    for (std::size_t i = 0; i < idxB_in.size(); ++i) packedB[i] = static_cast<std::uint8_t>(secret_vec[idxB_in[i]] - splitA);

    auto sharesA_packed = ss::getShares<std::uint8_t>(packedA, k, n, fieldA, kn);
    auto sharesB_packed = ss::getShares<std::uint8_t>(packedB, k, n, fieldB, kn);
    warnShareCount(sharesA_packed.size(), n);
    warnShareCount(sharesB_packed.size(), n);

    std::vector<std::vector<std::uint8_t>> shares(n, std::vector<std::uint8_t>(secret_vec.size(), 0));
    for (unsigned s = 0; s < n; ++s) {
        for (std::size_t j = 0; j < idxA_in.size(); ++j) shares[s][idxA_in[j]] = sharesA_packed[s][j];
        for (std::size_t j = 0; j < idxB_in.size(); ++j) shares[s][idxB_in[j]] = static_cast<std::uint8_t>(sharesB_packed[s][j] + splitA);
    }
    return shares;
}

template <typename FieldType>
static std::vector<std::uint8_t> reconstructAndSave(
    const std::vector<std::vector<std::uint8_t>>& orig_shares,
    const FieldType& field,
    const std::vector<unsigned>& indices,
    unsigned rec_k,
    unsigned kn,
    unsigned order,
    FieldBehavior behavior,
    const std::string& out_path,
    int width, int height,
    std::size_t secret_size
) {
    auto [selectedShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, orig_shares);
    auto rec = ss::reconstructFromShares<std::uint8_t>(selectedShares, selectedXs, rec_k, field, kn, static_cast<unsigned>(secret_size));

    if (behavior == FieldBehavior::STRETCH) {
        rec = ss::stretchHistogram(rec, 0, static_cast<std::uint8_t>(order - 1), 0, 255);
    }

    ss::saveGrayscalePNG(out_path, rec, width, height);
    return rec;
}

template <typename FieldType>
static std::vector<std::uint8_t> reconstructAndSaveKeepPlain(
    const std::vector<std::vector<std::uint8_t>>& shares_full,
    const FieldType& field,
    const std::vector<std::size_t>& in_idxs,
    const std::vector<unsigned>& indices,
    unsigned rec_k,
    unsigned kn,
    const std::vector<std::uint8_t>& secret_vec,
    const std::string& out_path,
    int width, int height
) {
    auto [selectedSharesFull, selectedXs] = ss::selectSharesAndEvalPoints(indices, shares_full);

    std::vector<std::vector<std::uint8_t>> selectedPacked(selectedSharesFull.size());
    for (std::size_t si = 0; si < selectedSharesFull.size(); ++si) {
        selectedPacked[si].reserve(in_idxs.size());
        for (std::size_t pos : in_idxs) selectedPacked[si].push_back(selectedSharesFull[si][pos]);
    }

    auto rec_packed = ss::reconstructFromShares<std::uint8_t>(selectedPacked, selectedXs, rec_k, field, kn, static_cast<unsigned>(in_idxs.size()));

    std::vector<std::uint8_t> rec_full(secret_vec.size(), 0);
    for (std::size_t j = 0; j < in_idxs.size() && j < rec_packed.size(); ++j) rec_full[in_idxs[j]] = rec_packed[j];

    std::vector<char> is_in(secret_vec.size(), 0);
    for (std::size_t p : in_idxs) is_in[p] = 1;

    for (std::size_t i = 0; i < secret_vec.size(); ++i) {
        if (!is_in[i]) {
            unsigned sum = 0;
            for (std::size_t si = 0; si < selectedSharesFull.size(); ++si) sum += selectedSharesFull[si][i];
            rec_full[i] = static_cast<std::uint8_t>((sum + selectedSharesFull.size()/2) / selectedSharesFull.size());
        }
    }

    ss::saveGrayscalePNG(out_path, rec_full, width, height);
    return rec_full;
}

static std::vector<std::uint8_t> reconstructAndSaveSplit(
    const std::vector<std::vector<std::uint8_t>>& shares_full,
    unsigned splitA, unsigned splitB,
    const std::vector<std::size_t>& idxA_in,
    const std::vector<std::size_t>& idxB_in,
    const std::vector<unsigned>& indices,
    unsigned rec_k,
    unsigned kn,
    const std::vector<std::uint8_t>& secret_vec,
    const std::string& out_path,
    int width, int height
) {
    auto [selectedFullShares, selectedXs] = ss::selectSharesAndEvalPoints(indices, shares_full);

    std::vector<std::vector<std::uint8_t>> selectedPackedA(selectedFullShares.size()), selectedPackedB(selectedFullShares.size());
    for (std::size_t si = 0; si < selectedFullShares.size(); ++si) {
        selectedPackedA[si].reserve(idxA_in.size());
        for (std::size_t pos : idxA_in) selectedPackedA[si].push_back(selectedFullShares[si][pos]);

        selectedPackedB[si].reserve(idxB_in.size());
        for (std::size_t pos : idxB_in) selectedPackedB[si].push_back(static_cast<std::uint8_t>(selectedFullShares[si][pos] - splitA));
    }

    ss::PrimeField<std::uint8_t> fieldA(splitA);
    ss::PrimeField<std::uint8_t> fieldB(splitB);

    auto rec_packedA = ss::reconstructFromShares<std::uint8_t>(selectedPackedA, selectedXs, rec_k, fieldA, kn, static_cast<unsigned>(idxA_in.size()));
    auto rec_packedB = ss::reconstructFromShares<std::uint8_t>(selectedPackedB, selectedXs, rec_k, fieldB, kn, static_cast<unsigned>(idxB_in.size()));

    std::vector<std::uint8_t> rec_full(secret_vec.size(), 0);
    for (std::size_t j = 0; j < idxA_in.size() && j < rec_packedA.size(); ++j) rec_full[idxA_in[j]] = rec_packedA[j];
    for (std::size_t j = 0; j < idxB_in.size() && j < rec_packedB.size(); ++j) rec_full[idxB_in[j]] = static_cast<std::uint8_t>(rec_packedB[j] + splitA);

    ss::saveGrayscalePNG(out_path, rec_full, width, height);
    return rec_full;
}

template <typename FieldType>
void runScheme(
    const std::string &label,
    FieldType field,
    unsigned n,
    unsigned k,
    unsigned kn,
    const std::vector<std::uint8_t> &secret_vec,
    std::size_t secret_size,
    int width,
    int height,
    const std::string &result_prefix,
    bool do_failed_experiment = false,
    bool divide_share_height_by_k = false,
    FieldBehavior behavior = FieldBehavior::KEEP_PLAIN
) {
    std::cout << "\n--- " << label << " (n=" << n << ", k=" << k << ", kn=" << kn << ") ---\n";

    fs::create_directories(result_prefix + "/shadows");

    const unsigned order = static_cast<unsigned>(field.getOrder());

    if (behavior == FieldBehavior::KEEP_PLAIN) {
        auto in_idxs = indicesWhere(secret_vec, [&](std::uint8_t v){ return v < order; });
        auto oob_idxs = indicesWhere(secret_vec, [&](std::uint8_t v){ return v >= order; });

        auto shares = createSharesKeepPlain(field, secret_vec, n, k, kn, in_idxs);

        for (unsigned s = 0; s < n; ++s) {
            for (std::size_t pos : oob_idxs) shares[s][pos] = secret_vec[pos];
        }

        saveSharesImages(result_prefix, shares, width, height, divide_share_height_by_k, k);

        auto idxA = makeIndicesA(k);
        auto idxB = makeIndicesB(static_cast<unsigned>(shares.size()), k);

        std::string recA_path = result_prefix + "/reconstruction1.png";
        auto recA = reconstructAndSaveKeepPlain(shares, field, in_idxs, idxA, k, kn, secret_vec, recA_path, width, height);
        printMetrics(recA_path, secret_vec, recA);

        std::string recB_path = result_prefix + "/reconstruction2.png";
        auto recB = reconstructAndSaveKeepPlain(shares, field, in_idxs, idxB, k, kn, secret_vec, recB_path, width, height);
        printMetrics(recB_path, secret_vec, recB);

        if (do_failed_experiment) {
            unsigned k_fail = (k > 1 ? k - 1 : 1);
            std::vector<unsigned> idxF;
            for (unsigned i = 0; i < k_fail && i < shares.size(); ++i) idxF.push_back(i);
            std::string recF_path = result_prefix + "/reconstruction_failed.png";
            auto recF = reconstructAndSaveKeepPlain(shares, field, in_idxs, idxF, k_fail, kn, secret_vec, recF_path, width, height);
            printMetrics(recF_path, secret_vec, recF);
        }

        if (kn == 1) {
            std::string share1_path = result_prefix + "/shadows/shadow_1.png";
            printMetrics(share1_path, secret_vec, shares[0]);
        }

        return;
    }
    else if (behavior == FieldBehavior::SPLIT) {
        const unsigned splitA = 107, splitB = 149;

        // Partition indices into two groups
        auto idxA_in = indicesWhere(secret_vec, [&](std::uint8_t v){ return v < splitA; });
        auto idxB_in = indicesWhere(secret_vec, [&](std::uint8_t v){ return v >= splitA; });

        auto shares = createSharesSplit(splitA, splitB, secret_vec, n, k, kn, idxA_in, idxB_in);

        saveSharesImages(result_prefix, shares, width, height, divide_share_height_by_k, k);

        auto idxA_sel = makeIndicesA(k);
        auto idxB_sel = makeIndicesB(static_cast<unsigned>(shares.size()), k);

        std::string recA_path = result_prefix + "/reconstruction1.png";
        auto recA = reconstructAndSaveSplit(shares, splitA, splitB, idxA_in, idxB_in, idxA_sel, k, kn, secret_vec, recA_path, width, height);
        printMetrics(recA_path, secret_vec, recA);

        std::string recB_path = result_prefix + "/reconstruction2.png";
        auto recB = reconstructAndSaveSplit(shares, splitA, splitB, idxA_in, idxB_in, idxB_sel, k, kn, secret_vec, recB_path, width, height);
        printMetrics(recB_path, secret_vec, recB);

        if (do_failed_experiment) {
            unsigned k_fail = (k > 1 ? k - 1 : 1);
            std::vector<unsigned> idxF;
            for (unsigned i = 0; i < k_fail && i < shares.size(); ++i) idxF.push_back(i);
            std::string recF_path = result_prefix + "/reconstruction_failed.png";
            auto recF = reconstructAndSaveSplit(shares, splitA, splitB, idxA_in, idxB_in, idxF, k_fail, kn, secret_vec, recF_path, width, height);
            printMetrics(recF_path, secret_vec, recF);
        }

        if (kn == 1) {
            std::string share1_path = result_prefix + "/shadows/shadow_1.png";
            printMetrics(share1_path, secret_vec, shares[0]);
        }

        return;
    }

    auto field_vec = preprocessFieldVec(secret_vec, field, behavior);
    auto orig_shares = createShares(field_vec, k, n, field, kn);
    saveSharesImages(result_prefix, orig_shares, width, height, divide_share_height_by_k, k);

    auto idxA = makeIndicesA(k);
    auto idxB = makeIndicesB(static_cast<unsigned>(orig_shares.size()), k);

    std::string recA_path = result_prefix + "/reconstruction1.png";
    auto recA = reconstructAndSave(orig_shares, field, idxA, k, kn, order, behavior, recA_path, width, height, secret_size);
    printMetrics(recA_path, secret_vec, recA);

    std::string recB_path = result_prefix + "/reconstruction2.png";
    auto recB = reconstructAndSave(orig_shares, field, idxB, k, kn, order, behavior, recB_path, width, height, secret_size);
    printMetrics(recB_path, secret_vec, recB);

    if (do_failed_experiment) {
        unsigned k_fail = (k > 1 ? k - 1 : 1);
        std::vector<unsigned> idxF;
        for (unsigned i = 0; i < k_fail && i < orig_shares.size(); ++i) idxF.push_back(i);
        std::string recF_path = result_prefix + "/reconstruction_failed.png";
        auto recF = reconstructAndSave(orig_shares, field, idxF, k_fail, kn, order, behavior, recF_path, width, height, secret_size);
        printMetrics(recF_path, secret_vec, recF);
    }

    if (kn == 1) {
        std::string share1_path = result_prefix + "/shadows/shadow_1.png";
        const auto &share1 = orig_shares[0];
        printMetrics(share1_path, secret_vec, share1);
    }
}

int main(int argc, char* argv[]) {
    const std::string input_path = "images/stop.png";
    auto [secret_pixels, width, height] = ss::readGrayscale(input_path);
    std::size_t secret_size = secret_pixels.size();

    const unsigned n = 6;
    const unsigned k = 4;

    FieldBehavior behavior;

    if (argc < 2) {
        std::cout << "0: CLAMP\n1: STRETCH\n2: KEEP_PLAIN\n3: SPLIT\n";
        return 1;
    } else {
        int option = std::stoi(argv[1]);
        switch (option) {
            case 0: behavior = FieldBehavior::CLAMP; break;
            case 1: behavior = FieldBehavior::STRETCH; break;
            case 2: behavior = FieldBehavior::KEEP_PLAIN; break;
            case 3: behavior = FieldBehavior::SPLIT; break;
            default:
                std::cerr << "Invalid option. Use 0, 1, 2, or 3.\n";
                return 1;
        }
    }

    // 1) Shamir GF(251)
    runScheme("Shamir GF(251)", ss::PrimeField<std::uint8_t>(251), n, k, 1, secret_pixels, secret_size, width, height, "results/shamir251", true, false, behavior);

    if (behavior == FieldBehavior::SPLIT) {
        return 0;
    }

    // 2) Shamir GF(2^8)
    runScheme("Shamir GF(2^8)", ss::BinaryField<std::uint8_t>(8), n, k, 1, secret_pixels, secret_size, width, height, "results/shamir256", true, false, behavior);

    if (behavior == FieldBehavior::KEEP_PLAIN) {
        return 0;
    }

    // 3) Thien-Lin GF(251)
    runScheme("Thien-Lin GF(251)", ss::PrimeField<std::uint8_t>(251), n, k, 4, secret_pixels, secret_size, width, height, "results/thienlin251", false, true, behavior);

    // 4) Thien-Lin GF(2^8)
    runScheme("Thien-Lin GF(2^8)", ss::BinaryField<std::uint8_t>(8), n, k, 4, secret_pixels, secret_size, width, height, "results/thienlin256", false, true, behavior);

    return 0;
}
