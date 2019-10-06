#pragma once

#include "trit_array_trie.hpp"
#include "vcode_array.hpp"

template <class Point, uint32_t Length>
class MultiTrieIndex {
  public:
    using point_type = Point;
    using vcode_trait_type = VCodeTraits<Length>;
    using vint_type = typename vcode_trait_type::vint_type;
    using hasher_type = FrechetHasher<Point, Length>;
    using database_type = VCodeArray<Point, Length>;
    using trie_type = TritArrayTrie;

  public:
    MultiTrieIndex() = default;

    MultiTrieIndex(const std::vector<Curve<Point>>& curves, const Configure& cfg)
        : m_indexes(cfg.buckets), m_bucket_begs(cfg.buckets + 1) {
        const size_t n_curves = curves.size();
        std::vector<uint8_t> codes_buf(n_curves * Length);
        std::vector<const uint8_t*> codes(n_curves);
        {
            size_t i = 0;
            m_database = std::make_unique<database_type>(curves, cfg, [&](const uint8_t* code) {
                const size_t beg = i * Length;
                std::copy(code, code + Length, &codes_buf[beg]);
                codes[i++] = &codes_buf[beg];
            });
        }

        uint32_t bkt_beg = 0;
        for (uint32_t b = 0; b < cfg.buckets; ++b) {
            m_bucket_begs[b] = bkt_beg;
            bkt_beg += (b + Length) / cfg.buckets;
        }
        m_bucket_begs[cfg.buckets] = bkt_beg;

        tfm::printfln("[MultiTrieIndex construction]");
        for (uint32_t b = 0; b < cfg.buckets; ++b) {
            uint32_t sub_length = m_bucket_begs[b + 1] - m_bucket_begs[b];
            tfm::printfln("# bucket %d (L=%d)", b, sub_length);
            m_indexes[b] = std::make_unique<trie_type>(codes, sub_length, cfg.reduction_threshold);
            for (size_t i = 0; i < n_curves; ++i) {
                codes[i] += sub_length;
            }
            tfm::printfln("# - %d nodes; %f values per leaf",  //
                          m_indexes[b]->get_num_nodes(), m_indexes[b]->get_num_values_per_leaf());
        }

        m_candidates.reserve(1U << 16);  // for range_search
    }

    ~MultiTrieIndex() = default;

    size_t range_search(const Curve<Point>& query_curve, uint32_t hamming_range, std::function<void(size_t)> fn) {
        m_candidates.clear();
        size_t num_candidates = 0;

        const uint8_t* qcode = m_database->get_hasher().hash_code8(query_curve.points);
        const vint_type* qvcode = vcode_trait_type::to_vints(qcode, DOMAIN_BITS);

        const int32_t buckets = static_cast<int32_t>(m_indexes.size());
        const int32_t gph_range = int32_t(hamming_range) - buckets + 1;

        for (int32_t b = 0; b < buckets; ++b) {
            if (gph_range + b < 0) {
                continue;
            }
            const int32_t sub_range = (gph_range + b) / buckets;
            const uint8_t* sub_qcode = qcode + m_bucket_begs[b];
            m_indexes[b]->range_search(sub_qcode, sub_range, [&](uint32_t id) { m_candidates.push_back(id); });
        }

        std::sort(m_candidates.begin(), m_candidates.end());

        uint32_t prev_id = uint32_t(-1);
        for (size_t i = 0; i < m_candidates.size(); ++i) {
            const uint32_t cand_id = m_candidates[i];
            if (prev_id == cand_id) {
                continue;
            }

            uint32_t hamming_dist =
                vcode_trait_type::get_hamming_distance(qvcode, m_database->get(cand_id), DOMAIN_BITS);
            if (hamming_dist <= hamming_range) {
                fn(cand_id);
            }

            prev_id = cand_id;
            ++num_candidates;
        }

        return num_candidates;
    }

    // in bytes
    size_t get_memory_usage() const {
        return get_index_memory_usage() + get_database_memory_usage();
    }
    size_t get_index_memory_usage() const {
        size_t memory = 0;
        for (size_t b = 0; b < m_indexes.size(); ++b) {
            memory += m_indexes[b]->get_memory_usage();
        }
        memory += m_bucket_begs.size() * sizeof(uint32_t);
        return memory;
    }
    size_t get_database_memory_usage() const {
        return m_database->get_memory_usage();
    }

    size_t get_ITLB() const {
        return get_index_ITLB() + get_database_memory_usage();
    }
    size_t get_index_ITLB() const {
        size_t memory = 0;
        for (size_t b = 0; b < m_indexes.size(); ++b) {
            memory += m_indexes[b]->get_ITLB();
        }
        memory += m_bucket_begs.size() * sizeof(uint32_t);
        return memory;
    }

    size_t size() const {
        return m_database->size();
    }
    const VCodeArray<Point, Length>& get_database() const {
        return *m_database.get();
    }

    size_t get_num_nodes() const {
        size_t num = 0;
        for (size_t b = 0; b < m_indexes.size(); ++b) {
            num += m_indexes[b]->get_num_nodes();
        }
        return num;
    }
    size_t get_num_inner_nodes() const {
        size_t num = 0;
        for (size_t b = 0; b < m_indexes.size(); ++b) {
            num += m_indexes[b]->get_num_inner_nodes();
        }
        return num;
    }
    float get_num_values_per_leaf() const {
        float sum = 0.0;
        for (size_t b = 0; b < m_indexes.size(); ++b) {
            sum += m_indexes[b]->get_num_values_per_leaf();
        }
        return sum / m_indexes.size();
    }

  private:
    std::unique_ptr<database_type> m_database;
    std::vector<std::unique_ptr<trie_type>> m_indexes;
    std::vector<uint32_t> m_bucket_begs;

    // for seach
    std::vector<uint32_t> m_candidates;
};