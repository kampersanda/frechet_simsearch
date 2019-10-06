#pragma once

#include "basic.hpp"

struct BitTools {
    static const uint8_t POPCNT_TABLE[256];

    static uint32_t popcnt(uint8_t x) {
        return POPCNT_TABLE[x];
    }
    static uint32_t popcnt(uint16_t x) {
        return POPCNT_TABLE[x & UINT8_MAX] + POPCNT_TABLE[x >> 8];
    }
    static uint32_t popcnt(uint32_t x) {
#ifdef __SSE4_2__
        return static_cast<uint32_t>(__builtin_popcount(x));
#else
        x = (x & 0x55555555U) + ((x & 0xAAAAAAAAU) >> 1);
        x = (x & 0x33333333U) + ((x & 0xCCCCCCCCU) >> 2);
        x = (x & 0x0F0F0F0FU) + ((x & 0xF0F0F0F0U) >> 4);
        x *= 0x01010101U;
        return x >> 24;
#endif
    }
    static uint32_t popcnt(uint64_t x) {
#ifdef __SSE4_2__
        return static_cast<uint32_t>(__builtin_popcountll(x));
#else
        x = (x & 0x5555555555555555ULL) + ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
        x = (x & 0x3333333333333333ULL) + ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
        x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
        x *= 0x0101010101010101ULL;
        return x >> 56;
#endif
    }

    static void set_bit(uint32_t& x, uint32_t i, bool bit = true) {
        assert(i < 32);
        if (bit) {
            x |= (1U << i);
        } else {
            x &= ~(1U << i);
        }
    }
    static void set_bit(uint64_t& x, uint32_t i, bool bit = true) {
        assert(i < 64);
        if (bit) {
            x |= (1ULL << i);
        } else {
            x &= ~(1ULL << i);
        }
    }

    static bool get_bit(uint32_t x, uint32_t i) {
        assert(i < 64);
        return (x & (1U << i)) != 0;
    }
    static bool get_bit(uint64_t x, uint32_t i) {
        assert(i < 64);
        return (x & (1ULL << i)) != 0;
    }
};

const uint8_t BitTools::POPCNT_TABLE[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2,
    3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3,
    3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
    6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4,
    3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4,
    5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6,
    6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
