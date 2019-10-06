#pragma once

#include <queue>

#include <sdsl/bit_vectors.hpp>

#include <succinctrits/rank_support.hpp>

#include "basic.hpp"

class TritArrayTrie {
  public:
    TritArrayTrie() = default;

    TritArrayTrie(std::vector<const uint8_t*>& codes, const uint32_t length, const uint32_t reduction_threshold = 1)
        : m_length(length), m_reduction_threshold(reduction_threshold) {
        std::vector<kv_type> kvs(codes.size());
        for (uint32_t i = 0; i < codes.size(); ++i) {
            kvs[i] = kv_type{codes[i], i};
        }
        std::sort(kvs.begin(), kvs.end(), [&](const kv_type& x, const kv_type& y) {
            return std::lexicographical_compare(x.code, x.code + m_length, y.code, y.code + m_length);
        });

        build_trie(kvs);
    }

    void find(const uint8_t* code, std::function<void(uint32_t)> fn) const {
        exact_search(code, 0, fn);
    }

    void range_search(const uint8_t* code, uint32_t hamming_range, std::function<void(uint32_t)> fn) const {
        range_search_recur(code, 0, hamming_range, fn);
    }

    size_t get_memory_usage() const {
        size_t memory = 0;
        memory += m_nodes.size_in_bytes();
        memory += m_nodes_r1.size_in_bytes();
        memory += m_nodes_r2.size_in_bytes();
        memory += sdsl::size_in_bytes(m_value_begs);
        memory += sdsl::size_in_bytes(m_value_begs_s1);
        memory += sdsl::size_in_bytes(m_values);
        memory += sizeof(m_length);
        memory += sizeof(m_reduction_threshold);
        memory += sizeof(m_num_inner_nodes);
        memory += sizeof(m_num_leaf_nodes);
        return memory;
    }

    size_t get_ITLB() const {
        size_t memory = 0;
        memory += size_t(9.44 * get_num_nodes() / 8ULL + 1ULL);
        memory += sdsl::size_in_bytes(m_value_begs);
        memory += sdsl::size_in_bytes(m_value_begs_s1);
        memory += sdsl::size_in_bytes(m_values);
        memory += sizeof(m_length);
        memory += sizeof(m_reduction_threshold);
        memory += sizeof(m_num_inner_nodes);
        memory += sizeof(m_num_leaf_nodes);
        return memory;
    }

    void debug_print(std::ostream& os) const {
        os << "-- TritArrayTrie --\n";

        std::vector<std::vector<uint32_t>> Ps;

        os << "NODES : ";
        for (uint32_t i = 0; i < m_nodes.get_num_trits(); i += MAX_FANOUT) {
            for (uint32_t j = i; j < i + MAX_FANOUT; ++j) {
                if (m_nodes[j] != 0) {
                    os << uint8_t(j % MAX_FANOUT) << "(" << int(m_nodes[j]) << ") ";
                }
            }
            os << "| ";
        }
        os << '\n';

        os << "VALUES: ";
        for (uint32_t i = 0; i < m_values.size(); ++i) {
            os << m_values[i] << " ";
            if (m_value_begs[i + 1]) {
                os << "| ";
            }
        }
        os << '\n';
    }

    uint32_t get_num_nodes() const {
        return m_num_inner_nodes + m_num_leaf_nodes;
    }
    uint32_t get_num_inner_nodes() const {
        return m_num_inner_nodes;
    }
    float get_num_values_per_leaf() const {
        return float(m_values.size()) / m_num_leaf_nodes;
    }

  private:
    static constexpr uint32_t MAX_FANOUT = 256;

    succinctrits::trit_vector m_nodes;  // L and C = 1, L and !C = 2
    succinctrits::rank_support<1> m_nodes_r1;
    succinctrits::rank_support<2> m_nodes_r2;
    sdsl::bit_vector m_value_begs;
    sdsl::bit_vector::select_1_type m_value_begs_s1;
    std::vector<uint32_t> m_values;
    uint32_t m_length = 0;
    uint32_t m_reduction_threshold = 0;
    uint32_t m_num_inner_nodes = 0;
    uint32_t m_num_leaf_nodes = 0;

    struct kv_type {
        const uint8_t* code;
        uint32_t id;
    };
    struct node_type {
        uint32_t beg;
        uint32_t end;
        uint32_t lev;
        uint32_t pos;
    };

    void build_trie(const std::vector<kv_type>& kvs) {
        std::queue<node_type> q;
        q.push(node_type{0, static_cast<uint32_t>(kvs.size()), 0, UINT32_MAX});

        std::vector<uint8_t> nodes;
        nodes.reserve(1U << 10);

        std::vector<bool> value_begs;
        std::vector<uint32_t> values;
        value_begs.reserve(kvs.size() + 1);
        values.reserve(kvs.size());

        m_num_inner_nodes = 0;
        m_num_leaf_nodes = 0;

        while (!q.empty()) {
            const uint32_t beg = q.front().beg;
            const uint32_t end = q.front().end;
            const uint32_t lev = q.front().lev;
            const uint32_t pos = q.front().pos;
            q.pop();

            const uint32_t num_values = end - beg;

            if (num_values <= m_reduction_threshold || lev == m_length) {
                assert(nodes[pos]);
                nodes[pos] = 2;  // does not have child

                for (uint32_t i = beg; i < end; ++i) {
                    value_begs.push_back(false);
                    values.push_back(kvs[i].id);
                }

                value_begs[value_begs.size() - (end - beg)] = true;
                ++m_num_leaf_nodes;

                continue;
            }

            if (nodes.size() == nodes.capacity()) {
                nodes.reserve(nodes.size() * 2);
            }
            nodes.resize(nodes.size() + MAX_FANOUT);

            const uint32_t child_spos = m_num_inner_nodes * MAX_FANOUT;
            assert(child_spos + MAX_FANOUT == nodes.size());

            uint8_t prev_c = kvs[beg].code[lev];
            uint32_t prev_beg = beg;

            for (uint32_t i = beg + 1; i < end; ++i) {
                const uint8_t curr_c = kvs[i].code[lev];
                assert(prev_c <= curr_c);
                if (prev_c != curr_c) {
                    const uint32_t child_pos = child_spos + prev_c;
                    nodes[child_pos] = 1;
                    q.push(node_type{prev_beg, i, lev + 1, child_pos});
                    prev_c = curr_c;
                    prev_beg = i;
                }
            }
            const uint32_t child_pos = child_spos + prev_c;
            nodes[child_pos] = 1;
            q.push(node_type{prev_beg, end, lev + 1, child_pos});

            ++m_num_inner_nodes;
        }

        value_begs.push_back(true);
        assert(values.size() == values.capacity());

        // Release
        m_nodes.build(nodes.begin(), nodes.size());
        m_nodes_r1.build(&m_nodes);
        m_nodes_r2.build(&m_nodes);

        m_value_begs = sdsl::bit_vector(value_begs.size());
        std::copy(value_begs.begin(), value_begs.end(), m_value_begs.begin());
        sdsl::util::init_support(m_value_begs_s1, &m_value_begs);

        m_values = std::move(values);
    }

    void exact_search(const uint8_t* code, uint32_t node_id, std::function<void(uint32_t)> fn) const {
        while (true) {
            const uint32_t child_pos = node_id * MAX_FANOUT + *code;
            const uint8_t child_trit = m_nodes[child_pos];

            if (child_trit == 0) {
                return;
            }

            if (child_trit == 2) {  // has child?
                uint32_t value_pos = m_value_begs_s1(m_nodes_r2(child_pos) + 1);
                do {
                    fn(m_values[value_pos]);
                } while (!m_value_begs[++value_pos]);

                return;
            }

            node_id = m_nodes_r1(child_pos) + 1;
            ++code;
        }
    }

    void range_search_recur(const uint8_t* code, uint32_t node_id, uint32_t hamming_range,
                            std::function<void(uint32_t)> fn) const {
        if (hamming_range == 0) {
            exact_search(code, node_id, fn);
            return;
        }

        uint32_t child_id = UINT32_MAX;
        uint32_t value_pos = UINT32_MAX;

        const uint32_t child_beg = node_id * MAX_FANOUT;

        for (uint32_t c = 0; c < MAX_FANOUT; ++c) {
            const uint32_t child_pos = child_beg + c;
            const uint8_t child_trit = m_nodes[child_pos];
            if (child_trit == 0) {  // not found?
                continue;
            }

            if (child_trit == 2) {  // has child?
                if (value_pos == UINT32_MAX) {
                    value_pos = m_value_begs_s1(m_nodes_r2(child_pos) + 1);
                }
                do {
                    fn(m_values[value_pos]);
                } while (!m_value_begs[++value_pos]);
                continue;
            }

            uint32_t next_hamming_range = hamming_range;
            if (*code != static_cast<uint8_t>(c)) {
                --next_hamming_range;
            }

            if (child_id == UINT32_MAX) {
                child_id = m_nodes_r1(child_pos);
            }

            range_search_recur(code + 1, ++child_id, next_hamming_range, fn);
        }
    }
};
