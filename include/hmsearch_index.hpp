#pragma once

#include <hmsearch/hmsearch.hpp>

#include "vcode_array.hpp"

namespace frechet_simsearch {

template <class Point, uint32_t Length>
class HmsearchIndex {
  public:
    using point_type = Point;
    using vcode_trait_type = VCodeTraits<Length>;
    using vint_type = typename vcode_trait_type::vint_type;
    using hasher_type = FrechetHasher<Point, Length>;

  public:
    HmsearchIndex() = default;

    HmsearchIndex(const std::vector<Curve<Point>>& curves, const Configure& cfg)
        : m_size(curves.size()), m_hasher(cfg, max_curve_length(curves)) {
        std::vector<uint8_t> codes(curves.size() * Length);
        for (size_t i = 0; i < curves.size(); ++i) {
            const uint8_t* code = m_hasher.hash_code8(curves[i].points);
            std::copy(code, code + Length, codes.data() + (i * Length));
        }
        std::vector<const uint8_t*> keys(curves.size());
        for (size_t i = 0; i < curves.size(); ++i) {
            keys[i] = codes.data() + (i * Length);
        }
        m_index = std::make_unique<hmsearch::hm_index>();
        m_index->build(keys, Length, 1U << DOMAIN_BITS, cfg.buckets);
    }

    ~HmsearchIndex() = default;

    size_t range_search(const Curve<Point>& query_curve, uint32_t hamming_range, std::function<void(uint32_t)> fn) {
        const uint8_t* qcode = m_hasher.hash_code8(query_curve.points);
        return m_index->search(qcode, hamming_range, fn);
    }

    size_t size() const {
        return m_size;
    }
    size_t get_memory_usage() const {
        return sdsl::size_in_bytes(*m_index.get());
    }

  private:
    size_t m_size = 0;
    hasher_type m_hasher;
    std::unique_ptr<hmsearch::hm_index> m_index;
};

}  // namespace frechet_simsearch