#pragma once

#include "bit_tools.hpp"

namespace frechet_simsearch {

template <class VInt, class Int, uint32_t Length>
inline VInt to_vint(const Int* in, uint32_t j) {
    VInt v = VInt(0);
    for (uint32_t i = 0; i < Length; ++i) {
        VInt b = (in[i] >> j) & VInt(1);
        v |= (b << i);
    }
    return v;
}

template <uint32_t Length>
struct VCodeTraits;

template <>
struct VCodeTraits<32> {
    using vint_type = uint32_t;
    static constexpr uint32_t LENGTH = 32;

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits) {
        vint_type diff = 0;
        for (uint32_t j = 0; j < domain_bits; ++j) {
            diff |= (x[j] ^ y[j]);
        }
        return BitTools::popcnt(diff);
    }

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits, uint32_t range) {
        vint_type diff = 0;
        uint32_t dist = 0;
        for (uint32_t j = 0; j < domain_bits; ++j) {
            diff |= (x[j] ^ y[j]);
            dist = BitTools::popcnt(diff);
            if (dist > range) {
                return dist;
            }
        }
        return dist;
    }

    template <class Int>
    static const vint_type* to_vints(const Int* in, uint32_t domain_bits) {
        static vint_type out[32];
        for (uint32_t j = 0; j < domain_bits; ++j) {
            out[j] = to_vint<vint_type, Int, LENGTH>(in, j);
        }
        return out;
    }

    template <class Int>
    static void to_vints(const Int* in, vint_type* out, uint32_t domain_bits) {
        for (uint32_t j = 0; j < domain_bits; ++j) {
            out[j] = to_vint<vint_type, Int, LENGTH>(in, j);
        }
    }
};

template <>
struct VCodeTraits<64> {
    using vint_type = uint64_t;
    static constexpr uint32_t LENGTH = 64;

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits) {
        vint_type diff = 0;
        for (uint32_t j = 0; j < domain_bits; ++j) {
            diff |= (x[j] ^ y[j]);
        }
        return BitTools::popcnt(diff);
    }

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits, uint32_t range) {
        vint_type diff = 0;
        uint32_t dist = 0;
        for (uint32_t j = 0; j < domain_bits; ++j) {
            diff |= (x[j] ^ y[j]);
            dist = BitTools::popcnt(diff);
            if (dist > range) {
                return dist;
            }
        }
        return dist;
    }

    template <class Int>
    static const vint_type* to_vints(const Int* in, uint32_t domain_bits) {
        static vint_type out[32];
        for (uint32_t j = 0; j < domain_bits; ++j) {
            out[j] = to_vint<vint_type, Int, LENGTH>(in, j);
        }
        return out;
    }

    template <class Int>
    static void to_vints(const Int* in, vint_type* out, uint32_t domain_bits) {
        for (uint32_t j = 0; j < domain_bits; ++j) {
            out[j] = to_vint<vint_type, Int, LENGTH>(in, j);
        }
    }
};

template <>
struct VCodeTraits<128> {
    using vint_type = std::tuple<uint64_t, uint64_t>;
    static constexpr uint32_t LENGTH = 128;

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits) {
        vint_type diff = {0, 0};
        for (uint32_t j = 0; j < domain_bits; ++j) {
            std::get<0>(diff) |= (std::get<0>(x[j]) ^ std::get<0>(y[j]));
            std::get<1>(diff) |= (std::get<1>(x[j]) ^ std::get<1>(y[j]));
        }
        return BitTools::popcnt(std::get<0>(diff)) + BitTools::popcnt(std::get<1>(diff));
    }

    static uint32_t get_hamming_distance(const vint_type* x, const vint_type* y, uint32_t domain_bits, uint32_t range) {
        uint32_t dist = 0;
        vint_type diff = {0, 0};
        for (uint32_t j = 0; j < domain_bits; ++j) {
            std::get<0>(diff) |= (std::get<0>(x[j]) ^ std::get<0>(y[j]));
            std::get<1>(diff) |= (std::get<1>(x[j]) ^ std::get<1>(y[j]));
            dist = BitTools::popcnt(std::get<0>(diff)) + BitTools::popcnt(std::get<1>(diff));
            if (dist > range) {
                return dist;
            }
        }
        return dist;
    }

    template <class Int>
    static const vint_type* to_vints(const Int* in, uint32_t domain_bits) {
        static vint_type out[32];
        for (uint32_t j = 0; j < domain_bits; ++j) {
            std::get<0>(out[j]) = to_vint<uint64_t, Int, 64>(in, j);
            std::get<1>(out[j]) = to_vint<uint64_t, Int, 64>(in + 64, j);
        }
        return out;
    }

    template <class Int>
    static void to_vints(const Int* in, vint_type* out, uint32_t domain_bits) {
        for (uint32_t j = 0; j < domain_bits; ++j) {
            std::get<0>(out[j]) = to_vint<uint64_t, Int, 64>(in, j);
            std::get<1>(out[j]) = to_vint<uint64_t, Int, 64>(in + 64, j);
        }
    }
};

}  // namespace frechet_simsearch