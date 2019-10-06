#pragma once

#include <bitset>

#include "frechet_hasher.hpp"
#include "vcode_traits.hpp"

template <class Point, uint32_t Length>
class VCodeArray {
  public:
    using point_type = Point;
    using vcode_trait_type = VCodeTraits<Length>;
    using vint_type = typename vcode_trait_type::vint_type;
    using hasher_type = FrechetHasher<Point, Length>;

  public:
    VCodeArray() = default;

    VCodeArray(const std::vector<Curve<Point>>& curves, const Configure& cfg)
        : m_hasher(cfg, max_curve_length(curves)), m_vcodes(curves.size() * DOMAIN_BITS) {
        tfm::printfln("[VCodeArray construction]");
        ProgressPrinter p(curves.size());
        for (size_t i = 0; i < curves.size(); ++i) {
            const vint_type* vcode = m_hasher.hash_vcode(curves[i].points);
            std::copy(vcode, vcode + DOMAIN_BITS, m_vcodes.data() + (i * DOMAIN_BITS));
            p(i + 1);
        }
    }

    VCodeArray(const std::vector<Curve<Point>>& curves, const Configure& cfg, std::function<void(const uint8_t*)> fn)
        : m_hasher(cfg, max_curve_length(curves)), m_vcodes(curves.size() * DOMAIN_BITS) {
        tfm::printfln("[VCodeArray construction]");
        ProgressPrinter p(curves.size() - 1);
        for (size_t i = 0; i < curves.size(); ++i) {
            const uint8_t* code = m_hasher.hash_code8(curves[i].points);
            fn(code);
            const vint_type* vcode = vcode_trait_type::to_vints(code, DOMAIN_BITS);
            std::copy(vcode, vcode + DOMAIN_BITS, m_vcodes.data() + (i * DOMAIN_BITS));
            p(i);
        }
    }

    ~VCodeArray() = default;

    // by LinearSearch
    void range_search(const Curve<Point>& query_curve, uint32_t hamming_range, std::function<void(size_t)> fn) const {
        const vint_type* qvcode = m_hasher.hash_vcode(query_curve.points);
        range_search(qvcode, hamming_range, fn);
    }
    void range_search(const vint_type* qvcode, uint32_t hamming_range, std::function<void(size_t)> fn) const {
        const size_t n_codes = size();
        for (size_t i = 0; i < n_codes; ++i) {
            uint32_t hamming_dist = vcode_trait_type::get_hamming_distance(qvcode, get(i), DOMAIN_BITS);
            if (hamming_dist <= hamming_range) {
                fn(i);
            }
        }
    }

    const vint_type* operator[](size_t i) const {
        return get(i);
    }
    const vint_type* get(size_t i) const {
        return &m_vcodes[i * DOMAIN_BITS];
    }

    const hasher_type& get_hasher() const {
        return m_hasher;
    }

    size_t get_memory_usage() const {
        size_t memory = 0;
        memory += sizeof(vint_type) * m_vcodes.size();
        return memory;
    }
    size_t size() const {
        return m_vcodes.size() / DOMAIN_BITS;
    }

  private:
    hasher_type m_hasher;
    std::vector<vint_type> m_vcodes;
};