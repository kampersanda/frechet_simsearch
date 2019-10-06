#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include <sdsl/int_vector.hpp>

// #define HMSEARCH_DISABLE_VERT
#define HMSEARCH_PRINT_PROGRESS

#define HMSEARCH_CHECK_IF(cond, msg)                    \
    do {                                                \
        if (cond) {                                     \
            std::cerr << "ERROR: " << msg << std::endl; \
            exit(1);                                    \
        }                                               \
    } while (0)

namespace hmsearch {

using signature_t = std::vector<uint32_t>;

inline bool operator==(const signature_t& x, const signature_t& y) {
    return std::equal(x.begin(), x.end(), y.begin(), y.end());
}
inline bool operator!=(const signature_t& x, const signature_t& y) {
    return !(x == y);
}

struct sig_hash {
    // fnv1a
    size_t operator()(const signature_t& sig) const {
        static const size_t init = size_t((sizeof(size_t) == 8) ? 0xcbf29ce484222325 : 0x811c9dc5);
        static const size_t multiplier = size_t((sizeof(size_t) == 8) ? 0x100000001b3 : 0x1000193);
        size_t h = init;
        for (size_t i = 0; i < sig.size(); ++i) {
            h ^= sig[i];
            h *= multiplier;
        }
        return h;
    }
    static const sig_hash& get_instance() {
        static sig_hash hasher;
        return hasher;
    }
};

#ifdef HMSEARCH_PRINT_PROGRESS
class progress_printer {
  public:
    progress_printer(size_t max) {
        const size_t step = max / 10;
        for (size_t i = 1; i <= 9; ++i) {
            m_stations[i - 1] = i * step;
        }
        m_stations[9] = max;
        m_stations[10] = size_t(-1);
    }

    void operator()(size_t i) {
        while (m_stations[m_cursor] <= i) {
            std::cerr << " *" << std::flush;
            ++m_cursor;
        }
        if (m_cursor == 10) {
            std::cerr << " done!!" << std::endl;
        }
    }

  private:
    size_t m_stations[11];
    size_t m_cursor = 0;
};
#endif

// one-del-var
class odv_index {
  public:
    using size_type = uint64_t;  // for sdsl::serialize

  private:
    static constexpr float LOAD_FACTOR = 1.5;

    struct element_t {
        uint32_t sig_pos;
        uint32_t id_beg;
        uint32_t id_end;
    };
    std::vector<element_t> m_table;
    std::vector<uint32_t> m_ids;
    sdsl::int_vector<> m_signatures;
    uint32_t m_length = 0;
    uint32_t m_del_marker = 0;

  public:
    odv_index() = default;
    ~odv_index() = default;

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_table, out);
        written_bytes += sdsl::serialize(m_ids, out);
        written_bytes += sdsl::serialize(m_signatures, out);
        written_bytes += sdsl::serialize(m_length, out);
        written_bytes += sdsl::serialize(m_del_marker, out);
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_table, in);
        sdsl::load(m_ids, in);
        sdsl::load(m_signatures, in);
        sdsl::load(m_length, in);
        sdsl::load(m_del_marker, in);
    }

    template <class T>
    void build(const std::vector<const T*>& keys, uint32_t length, uint32_t alphabet_size) {
        static_assert(sizeof(T) <= 4, "");

        HMSEARCH_CHECK_IF(alphabet_size == UINT32_MAX, "alphabet_size is too large.");
        HMSEARCH_CHECK_IF(keys.size() * m_length > UINT32_MAX, "size of ids exceeds.");

        m_length = length;
        m_del_marker = alphabet_size;

        std::unordered_map<signature_t, std::vector<uint32_t>, sig_hash> signature_map;
        {
#ifdef HMSEARCH_PRINT_PROGRESS
            std::cerr << " #    - Making signatures... " << std::flush;
            progress_printer p(keys.size() - 1);
#endif
            signature_t sig(m_length);
            for (uint32_t i = 0; i < keys.size(); ++i) {
                for (uint32_t j = 0; j < m_length; ++j) {
                    HMSEARCH_CHECK_IF(keys[i][j] >= alphabet_size, "keys include a large character.");

                    make_signature(keys[i], j, sig);

                    auto it = signature_map.find(sig);
                    if (it != signature_map.end()) {
                        it->second.push_back(i);
                    } else {
                        signature_map.insert(std::make_pair(sig, std::vector<uint32_t>{i}));
                    }
                }
#ifdef HMSEARCH_PRINT_PROGRESS
                p(i);
#endif
            }
        }

        HMSEARCH_CHECK_IF(signature_map.size() > UINT32_MAX, "number of signatures exceeds UINT32_MAX.");

        const size_t table_size = static_cast<size_t>(signature_map.size() * LOAD_FACTOR);
        m_table.resize(table_size, element_t{UINT32_MAX, 0, 0});
        m_ids.reserve(keys.size() * m_length);
        m_signatures = sdsl::int_vector<>(signature_map.size() * m_length, 0, sdsl::bits::hi(alphabet_size) + 1);

#ifdef HMSEARCH_PRINT_PROGRESS
        std::cerr << " #    - Storing signatures..." << std::flush;
        uint32_t progress = 0;
        progress_printer p(signature_map.size() - 1);
#endif

        uint64_t sig_beg = 0;

        for (const auto& kv : signature_map) {
            const signature_t& sig = kv.first;
            const std::vector<uint32_t>& ids = kv.second;

            uint64_t pos = sig_hash::get_instance()(sig) % table_size;

            while (true) {
                if (m_table[pos].sig_pos == UINT32_MAX) {  // vacant?
                    m_table[pos].sig_pos = sig_beg / m_length;
                    std::copy(sig.begin(), sig.end(), m_signatures.begin() + sig_beg);
                    sig_beg += m_length;  // sig.size()

                    m_table[pos].id_beg = m_ids.size();
                    std::copy(ids.begin(), ids.end(), std::back_inserter(m_ids));
                    m_table[pos].id_end = m_ids.size();

                    break;
                }

                ++pos;
                if (pos == table_size) {
                    pos = 0;
                }
            }
#ifdef HMSEARCH_PRINT_PROGRESS
            p(progress++);
#endif
        }

        assert(sig_beg == m_signatures.size());
        assert(m_ids.size() == keys.size() * m_length);
    }

    template <class T>
    void search(const T* key, signature_t& sig, std::function<void(uint32_t)> fn) const {
        sig.resize(m_length);

        for (uint32_t j = 0; j < m_length; ++j) {
            make_signature(key, j, sig);
            uint64_t pos = sig_hash::get_instance()(sig) % m_table.size();

            while (true) {
                if (m_table[pos].sig_pos == UINT32_MAX) {  // vacant?
                    break;
                }

                const uint64_t sig_beg = m_table[pos].sig_pos * m_length;

                if (std::equal(sig.begin(), sig.end(), m_signatures.begin() + sig_beg)) {
                    for (uint32_t i = m_table[pos].id_beg; i < m_table[pos].id_end; ++i) {
                        fn(m_ids[i]);
                    }
                    break;
                }

                ++pos;
                if (pos == m_table.size()) {
                    pos = 0;
                }
            }
        }
    }

    template <class T>
    void make_signature(const T* key, uint32_t i, signature_t& out) const {
        assert(out.size() == m_length);
        std::copy(key, key + m_length, out.begin());
        out[i] = m_del_marker;
    }
};

class hm_index {
  public:
    using size_type = uint64_t;  // for sdsl::serialize

  private:
    std::vector<odv_index> m_odv_indexes;
    std::vector<uint32_t> m_bucket_begs;
    uint32_t m_length = 0;
    uint32_t m_alphabet_size = 0;
    uint32_t m_buckets = 0;
#ifdef HMSEARCH_DISABLE_VERT
    sdsl::int_vector<> m_keys;
#else
    sdsl::int_vector<> m_vertical_keys;
    uint32_t m_vertical_levels = 0;
#endif

  public:
    hm_index() = default;
    ~hm_index() = default;

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const {
        auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::serialize(m_odv_indexes, out);
        written_bytes += sdsl::serialize(m_bucket_begs, out);
        written_bytes += sdsl::serialize(m_length, out);
        written_bytes += sdsl::serialize(m_alphabet_size, out);
        written_bytes += sdsl::serialize(m_buckets, out);
#ifdef HMSEARCH_DISABLE_VERT
        written_bytes += sdsl::serialize(m_keys, out);
#else
        written_bytes += sdsl::serialize(m_vertical_keys, out);
        written_bytes += sdsl::serialize(m_vertical_levels, out);
#endif
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream& in) {
        sdsl::load(m_odv_indexes, in);
        sdsl::load(m_bucket_begs, in);
        sdsl::load(m_length, in);
        sdsl::load(m_alphabet_size, in);
        sdsl::load(m_buckets, in);
#ifdef HMSEARCH_DISABLE_VERT
        sdsl::load(m_keys, in);
#else
        sdsl::load(m_vertical_keys, in);
        sdsl::load(m_vertical_levels, in);
#endif
    }

    uint32_t get_length() const {
        return m_length;
    }
    uint32_t get_alphabet_size() const {
        return m_alphabet_size;
    }
    uint32_t get_buckets() const {
        return m_buckets;
    }
#ifndef HMSEARCH_DISABLE_VERT
    uint32_t get_vertical_levels() const {
        return m_vertical_levels;
    }
#endif

    static uint32_t get_proper_buckets(uint32_t range) {
        return (range + 3) / 2;
    }

    template <class T>
    void build(const std::vector<const T*>& keys, uint32_t length, uint32_t alphabet_size, uint32_t buckets) {
        HMSEARCH_CHECK_IF(length > 64, "length > 64 is not supported.");

#ifdef HMSEARCH_PRINT_PROGRESS
        std::cerr << " # [hm_index::build] buckets = " << buckets << std::endl;
#endif

        m_length = length;
        m_alphabet_size = alphabet_size;
        m_buckets = buckets;

        m_odv_indexes.resize(m_buckets);
        m_bucket_begs.resize(m_buckets + 1);

        uint32_t bucket_beg = 0;
        for (uint32_t b = 0; b < m_buckets; ++b) {
            m_bucket_begs[b] = bucket_beg;
            bucket_beg += (m_length + b) / m_buckets;
        }
        m_bucket_begs[m_buckets] = bucket_beg;

        std::vector<const T*> bucket_keys(keys.size());
        for (uint32_t b = 0; b < m_buckets; ++b) {
#ifdef HMSEARCH_PRINT_PROGRESS
            std::cerr << " #   - bucket_id = " << b << std::endl;
#endif
            for (size_t i = 0; i < keys.size(); ++i) {
                bucket_keys[i] = keys[i] + m_bucket_begs[b];
            }
            m_odv_indexes[b].build(bucket_keys, m_bucket_begs[b + 1] - m_bucket_begs[b], alphabet_size);
        }

#ifdef HMSEARCH_DISABLE_VERT
        m_keys = sdsl::int_vector<>(keys.size() * m_length, 0, sdsl::bits::hi(alphabet_size) + 1);
        for (size_t i = 0; i < keys.size(); ++i) {
            std::copy(keys[i], keys[i] + m_length, m_keys.begin() + (i * m_length));
        }
#else
        m_vertical_levels = sdsl::bits::hi(alphabet_size) + 1;
        m_vertical_keys = sdsl::int_vector<>(keys.size() * m_vertical_levels, 0, m_length);
        for (size_t i = 0; i < keys.size(); ++i) {
            const size_t beg = i * m_vertical_levels;
            for (uint32_t j = 0; j < m_vertical_levels; ++j) {
                m_vertical_keys[beg + j] = make_vertical_code(keys[i], m_length, j);
            }
        }
#endif
    }

    template <class T>
    uint64_t search(const T* query, uint32_t hamming_range, std::function<void(uint32_t)> fn) const {
        HMSEARCH_CHECK_IF(m_buckets != get_proper_buckets(hamming_range), "unsupported hamming range.");

        signature_t sig;
        std::unordered_map<uint32_t, uint32_t> match_map;
        std::unordered_map<uint32_t, std::vector<uint32_t>> cand_map;

        for (uint32_t b = 0; b < m_buckets; ++b) {
            const T* b_query = query + m_bucket_begs[b];
            const odv_index& odv_idx = m_odv_indexes[b];

            match_map.clear();

            odv_idx.search(b_query, sig, [&](uint32_t id) {
                auto it = match_map.find(id);
                if (it == match_map.end()) {
                    match_map.insert(std::make_pair(id, 1U));
                } else {
                    it->second += 1;
                }
            });

            for (const auto& kv : match_map) {
                if (kv.second > 2) {
                    auto it = cand_map.find(kv.first);
                    if (it != cand_map.end()) {
                        it->second.push_back(0);
                    } else {
                        cand_map.insert(std::make_pair(kv.first, std::vector<uint32_t>{0}));
                    }
                } else {
                    auto it = cand_map.find(kv.first);
                    if (it != cand_map.end()) {
                        it->second.push_back(1);
                    } else {
                        cand_map.insert(std::make_pair(kv.first, std::vector<uint32_t>{1}));
                    }
                }
            }
        }

#ifndef HMSEARCH_DISABLE_VERT
        std::vector<uint64_t> vertical_query(m_vertical_levels);
        for (uint32_t j = 0; j < m_vertical_levels; ++j) {
            vertical_query[j] = make_vertical_code(query, m_length, j);
        }
#endif

        uint64_t num_candidates = 0;

        for (const auto& kv : cand_map) {
            uint32_t cand_id = kv.first;
            const std::vector<uint32_t>& errors = kv.second;
            assert(errors.size() > 0);

            // enhanced filter
            bool filtered = false;

            if (hamming_range % 2 == 0) {
                if (errors.size() < 2) {  // has less than two number
                    if (errors[0] == 1) {
                        filtered = true;
                    }
                }
            } else {
                if (errors.size() < 3) {  // has less than three number
                    if (errors.size() == 1) {
                        filtered = true;
                    } else if (errors[0] == 1 && errors[1] == 1) {
                        filtered = true;
                    }
                }
            }

            // verification
            if (!filtered) {
                uint32_t hammina_dist = 0;
#ifdef HMSEARCH_DISABLE_VERT
                auto key = m_keys.begin() + (cand_id * m_length);
                for (uint32_t j = 0; j < m_length; ++j) {
                    if (query[j] != key[j]) {
                        ++hammina_dist;
                        if (hammina_dist > hamming_range) {
                            break;
                        }
                    }
                }
#else
                uint64_t cumdiff = 0;
                uint64_t beg = cand_id * m_vertical_levels;
                for (uint32_t j = 0; j < m_vertical_levels; ++j) {
                    uint64_t diff = m_vertical_keys[beg + j] ^ vertical_query[j];
                    cumdiff |= diff;
                    hammina_dist = sdsl::bits::cnt(cumdiff);
                    if (hammina_dist > hamming_range) {
                        break;
                    }
                }
#endif
                if (hammina_dist <= hamming_range) {
                    fn(cand_id);
                }

                ++num_candidates;
            }
        }

        return num_candidates;
    }

  private:
    template <class T>
    static uint64_t make_vertical_code(const T* key, uint32_t length, uint32_t level) {
        assert(length <= 64);

        uint64_t code = 0;
        for (uint32_t j = 0; j < length; ++j) {
            uint64_t bit = (key[j] >> level) & 1ULL;
            code |= (bit << j);
        }
        return code;
    }
};

}  // namespace hmsearch