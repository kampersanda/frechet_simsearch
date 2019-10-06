#pragma once

#include <FRESH/core/query.h>

#include "frechet_tools.hpp"
#include "vcode_traits.hpp"

template <typename Point, uint32_t Length>
class FrechetHasher {
  public:
    using vcode_trait_type = VCodeTraits<Length>;
    using vint_type = typename vcode_trait_type::vint_type;

  public:
    FrechetHasher() = default;

    FrechetHasher(const Configure& cfg, size_t max_curve_length) {
        Splitmix64 seeder(cfg.seed);
        Xorshift1024star rnd(seeder.next());

        m_functions.resize(Length);
        for (uint32_t i = 0; i < Length; ++i) {
            m_functions[i] = LinearFactorLSH<Point>(cfg.resolution, max_curve_length, rnd);
        }
    }

    const uint8_t* hash_code8(const std::vector<Point>& input) const {
        static uint8_t code[Length];
        for (uint32_t i = 0; i < Length; ++i) {
            code[i] = static_cast<uint8_t>(m_functions[i].hash32(input) & 0xFF);
        }
        return code;
    }
    const uint32_t* hash_code32(const std::vector<Point>& input) const {
        static uint32_t code[Length];
        for (uint32_t i = 0; i < Length; ++i) {
            code[i] = static_cast<uint32_t>(m_functions[i].hash32(input));
        }
        return code;
    }
    const vint_type* hash_vcode(const std::vector<Point>& input) const {
        return vcode_trait_type::to_vints(hash_code32(input), DOMAIN_BITS);
    }

  private:
    std::vector<LinearFactorLSH<Point>> m_functions;
};
