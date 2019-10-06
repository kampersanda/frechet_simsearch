/*
 * Copyright 2018 Matteo Ceccarello
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
#pragma once

#ifndef CURVEDIST_HASH_H
#define CURVEDIST_HASH_H

#include <boost/functional/hash.hpp>
#include "frechet.h"
#include "rand.h"
#include "threading.h"
#include "timer.h"
#include "types.h"

typedef std::pair<size_t, size_t> curve_to_idx_t;

extern std::vector<size_t> g_tl_tensored_multihash_scratch_left;
extern std::vector<size_t> g_tl_tensored_multihash_scratch_right;
#pragma omp threadprivate(g_tl_tensored_multihash_scratch_left, g_tl_tensored_multihash_scratch_right)

extern curve_to_idx_t g_tl_scratch_curve_to_idx;
#pragma omp threadprivate(g_tl_scratch_curve_to_idx)

extern std::vector<size_t> g_tl_hash_values;
#pragma omp threadprivate(g_tl_hash_values)

template <typename Point>
Point random_shift(double resolution, Xorshift1024star& rnd);

template <>
inline Point1D random_shift(double range, Xorshift1024star& rnd) {
    return Point1D(rnd.next_double(range));
}

template <>
inline Point2D random_shift(double range, Xorshift1024star& rnd) {
    return Point2D(rnd.next_double(range), rnd.next_double(range));
}

template <typename Point>
size_t num_random_coefficients(size_t num);
template <>
inline size_t num_random_coefficients<Point1D>(size_t num) {
    return num;
};
template <>
inline size_t num_random_coefficients<Point2D>(size_t num) {
    return 2 * num;
};

template <typename Point>
struct LinearFactorLSH {
  public:
    // Added by Shunsuke Kanda (https://github.com/kampersanda)
    LinearFactorLSH() = default;

    LinearFactorLSH(double resolution, size_t max_curve_len, Xorshift1024star& rnd)
        : m_grid_delta(resolution),
          m_shift(random_shift<Point>(resolution, rnd)),
          m_rand_coefficients(num_random_coefficients<Point>(max_curve_len)) {
        for (size_t i = 0; i < m_rand_coefficients.size(); i++) {
            m_rand_coefficients[i] = rnd.next();
        }
    }

    void hash(const std::vector<Point>& input_curve, std::vector<Point>& output_hash) const;

    size_t hash32(const std::vector<Point>& input_curve) const;

  private:
    double m_grid_delta;
    Point m_shift;
    std::vector<int64_t> m_rand_coefficients;
};

template <>
inline void LinearFactorLSH<Point2D>::hash(const std::vector<Point2D>& input_curve,
                                           std::vector<Point2D>& output_hash) const {
    size_t n = input_curve.size();

    double last_x = std::numeric_limits<double>::infinity(), last_y = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < n; i++) {
        double x = std::round((input_curve[i].x + m_shift.x) / m_grid_delta);
        double y = std::round((input_curve[i].y + m_shift.y) / m_grid_delta);
        if (x != last_x || y != last_y) {
            last_x = x;
            last_y = y;
            output_hash.emplace_back(x, y);
        }
    }
}

template <>
inline void LinearFactorLSH<Point1D>::hash(const std::vector<Point1D>& input_curve,
                                           std::vector<Point1D>& output_hash) const {
    size_t n = input_curve.size();

    double last_x = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < n; i++) {
        double x = std::round((input_curve[i].x + m_shift.x) / m_grid_delta);
        if (x != last_x) {
            last_x = x;
            output_hash.emplace_back(x);
        }
    }
}

template <>
inline size_t LinearFactorLSH<Point1D>::hash32(const std::vector<Point1D>& input_curve) const {
    size_t n = input_curve.size();

    int32_t last_x = std::numeric_limits<int32_t>::max();

    int64_t sum = 0;
    size_t coeff_idx = 0;

    for (size_t i = 0; i < n; i++) {
        // FIXME: Some case must be taken here with rounding.
        int32_t x = std::lround((input_curve[i].x + m_shift.x) / m_grid_delta);
        if (x != last_x) {
            last_x = x;
            sum += x * m_rand_coefficients[coeff_idx++];
        }
    }

    return sum >> 32;
}

template <>
inline size_t LinearFactorLSH<Point2D>::hash32(const std::vector<Point2D>& input_curve) const {
    size_t n = input_curve.size();

    int32_t last_x = std::numeric_limits<int32_t>::max();
    int32_t last_y = std::numeric_limits<int32_t>::max();

    uint64_t sum = 0;

    size_t coeff_idx = 0;

    for (size_t i = 0; i < n; i++) {
        // FIXME: Some case must be taken here with rounding.
        int32_t x = std::lround((input_curve[i].x + m_shift.x) / m_grid_delta);
        int32_t y = std::lround((input_curve[i].y + m_shift.y) / m_grid_delta);
        if (x != last_x || y != last_y) {
            last_x = x;
            last_y = y;
            sum += x * m_rand_coefficients[coeff_idx] + y * m_rand_coefficients[coeff_idx + 1];
            coeff_idx += 2;
        }
    }

    return sum >> 32;
}

template <typename Point>
class TensoredMultiHash {
  public:
    TensoredMultiHash(size_t k, size_t L, double resolution, size_t max_curve_len, Xorshift1024star& rnd)
        : m_left_functions(std::vector<std::vector<LinearFactorLSH<Point>>>()),
          m_right_functions(std::vector<std::vector<LinearFactorLSH<Point>>>()),
          m_k(k),
          m_repetitions(L),
          m_rand_coefficients(k) {
        const size_t sqrt_L = std::ceil(std::sqrt(L));
        const size_t k_left = std::ceil(k / 2.0);
        const size_t k_right = std::floor(k / 2.0);

        m_left_functions.reserve(sqrt_L);
        m_right_functions.reserve(sqrt_L);

        for (size_t i = 0; i < sqrt_L; i++) {
            std::vector<LinearFactorLSH<Point>> concats_left;
            for (size_t j = 0; j < k_left; j++) {
                concats_left.emplace_back(resolution, max_curve_len, rnd);
            }
            m_left_functions.push_back(concats_left);

            std::vector<LinearFactorLSH<Point>> concats_right;
            for (size_t j = 0; j < k_right; j++) {
                concats_right.emplace_back(resolution, max_curve_len, rnd);
            }
            m_right_functions.push_back(concats_right);
        }

        for (size_t i = 0; i < m_rand_coefficients.size(); i++) {
            m_rand_coefficients[i] = rnd.next();
        }
    }

    void hash(const std::vector<Point>& input_curve, std::vector<size_t>& output) {
        output.clear();
        output.reserve(m_repetitions);

        g_tl_tensored_multihash_scratch_left.clear();
        g_tl_tensored_multihash_scratch_right.clear();
        for (const auto fns_left : m_left_functions) {
            uint64_t sum = 0;
            size_t coeff_idx = 0;
            for (const auto fn : fns_left) {
                sum += fn.hash32(input_curve) * m_rand_coefficients[coeff_idx++];
            }
            size_t hash_value = sum >> 32;
            g_tl_tensored_multihash_scratch_left.push_back(hash_value);
        }
        for (const auto fns_right : m_right_functions) {
            uint64_t sum = 0;
            size_t coeff_idx = 0;
            for (const auto fn : fns_right) {
                sum += fn.hash32(input_curve) * m_rand_coefficients[coeff_idx++];
            }
            size_t hash_value = sum >> 32;
            g_tl_tensored_multihash_scratch_right.push_back(hash_value);
        }

        size_t rep_idx = 0;
        for (const size_t hash64_left : g_tl_tensored_multihash_scratch_left) {
            for (const size_t hash64_right : g_tl_tensored_multihash_scratch_right) {
                if (rep_idx >= m_repetitions) {
                    return;
                }
                int64_t hash64 = hash64_left * m_rand_coefficients[0] + hash64_right * m_rand_coefficients[1];
                size_t hash32 = hash64 >> 32;
                output.push_back(hash32);
                rep_idx++;
            }
        }
    }

  private:
    std::vector<std::vector<LinearFactorLSH<Point>>> m_left_functions;
    std::vector<std::vector<LinearFactorLSH<Point>>> m_right_functions;
    size_t m_k;
    size_t m_repetitions;
    std::vector<uint64_t> m_rand_coefficients;
};

template <typename Point>
class LSHTable {
  public:
    LSHTable(size_t num_curves, size_t max_curve_length)
        : m_hash_values(std::vector<curve_to_idx_t>(num_curves)),
          m_bucket_idx(),
          m_max_curve_length(max_curve_length),
          m_fixed(false) {}

    // void put_all(const std::vector<Curve<Point>> &dataset);

    inline void put_single(const size_t curve_id, const size_t hash64) {
        curve_to_idx_t& p = m_hash_values[curve_id];
        p.first = hash64;
        p.second = curve_id;
    }

    /* Fixes the content of the table, by sorting and indexing it */
    void fix_table() {
        START_TIMER(table_sorting);
        std::sort(m_hash_values.begin(), m_hash_values.end());
        STOP_TIMER(table_sorting);

        START_TIMER(bucket_indexing);
        size_t number_distinct = 1;
        size_t begin_idx = 0, end_idx = 0;
        size_t& last = m_hash_values[0].first;
        for (size_t i = 1; i <= m_hash_values.size(); i++) {
            if (i == m_hash_values.size() || last != m_hash_values[i].first) {
                number_distinct++;
                end_idx = i;
                m_bucket_idx[last] = std::make_pair(begin_idx, end_idx);
                begin_idx = i;
                if (i != m_hash_values.size()) {
                    last = m_hash_values.at(i).first;
                }
            }
        }
        STOP_TIMER(bucket_indexing);
        /* printf("There are %lu distinct hash values\n", number_distinct); */
        m_fixed = true;
    }

    /**
     * The collision counters refer to the curves passed to the `put` call, in
     * the same order.
     */
    void count_collisions(const size_t hash64, std::vector<size_t>& counters) {
        if (!m_fixed) {
            throw std::logic_error("The table has not been fixed");
        }
        std::pair<size_t, size_t> bucket_bounds = m_bucket_idx[hash64];
        for (size_t i = bucket_bounds.first; i < bucket_bounds.second; i++) {
            size_t idx = m_hash_values[i].second;
            counters.at(idx)++;
        }
    }

    void foreach_collision(const size_t hash64, std::vector<bool>& visited, std::function<void(size_t)> fn) {
        if (!m_fixed) {
            throw std::logic_error("The table has not been fixed");
        }
        std::pair<size_t, size_t> bucket_bounds = m_bucket_idx[hash64];
        for (size_t i = bucket_bounds.first; i < bucket_bounds.second; i++) {
            size_t idx = m_hash_values[i].second;
            if (!visited[idx]) {
                visited[idx] = true;
                fn(idx);
            }
        }
    }

  private:
    std::vector<curve_to_idx_t> m_hash_values;
    /// Stores begin/end indices (resp. inclusive and exclusive) of buckets in the (sorted) m_hash_values vector
    std::unordered_map<size_t, std::pair<size_t, size_t>> m_bucket_idx;
    size_t m_max_curve_length;
    bool m_fixed;
};

template <typename Point>
class MultiTable {
  public:
    MultiTable(size_t k, size_t L, double resolution, size_t num_curves, size_t max_curve_length, Xorshift1024star& rnd)
        : k(k),
          L(L),
          m_hash_functions(TensoredMultiHash<Point>(k, L, resolution, max_curve_length, rnd)),
          m_max_curve_length(max_curve_length) {
        m_tables.reserve(L);
        for (size_t i = 0; i < L; i++) {
            m_tables.emplace_back(num_curves, max_curve_length);
        }
    }

    void put_all(const std::vector<Curve<Point>>& dataset) {
        const size_t n = dataset.size();
        START_TIMER(table_hashing);
#pragma omp parallel for
        for (size_t i = 0; i < n; i++) {
            g_tl_hash_values.clear();
            m_hash_functions.hash(dataset[i].points, g_tl_hash_values);
            for (size_t rep_idx = 0; rep_idx < L; rep_idx++) {
                m_tables[rep_idx].put_single(i, g_tl_hash_values[rep_idx]);
            }
        }
        STOP_TIMER(table_hashing);

        // TODO: This can possibly be done in parallel
        for (auto& t : m_tables) {
            t.fix_table();
        }
    }

    void count_collisions(const std::vector<Point>& query_curve, std::vector<size_t>& counters) {
        std::fill(counters.begin(), counters.end(), 0);
        g_tl_hash_values.clear();
        m_hash_functions.hash(query_curve, g_tl_hash_values);
        for (size_t rep_idx = 0; rep_idx < L; rep_idx++) {
            m_tables[rep_idx].count_collisions(g_tl_hash_values[rep_idx], counters);
        }
    }

    void foreach_collision(const std::vector<Point>& query_curve, std::vector<bool>& visited,
                           std::function<void(size_t)> fn) {
        std::fill(visited.begin(), visited.end(), false);
        g_tl_hash_values.clear();
        m_hash_functions.hash(query_curve, g_tl_hash_values);
        for (size_t rep_idx = 0; rep_idx < L; rep_idx++) {
            m_tables[rep_idx].foreach_collision(g_tl_hash_values[rep_idx], visited, fn);
        }
    }

  private:
    /// The length of each hash (i.e. the number of concatenations of each hash)
    size_t k;
    /// The number of tables to use
    size_t L;
    /// The multiple hash functions (with repetitions) to use in to hash curves
    TensoredMultiHash<Point> m_hash_functions;
    /// The tables
    std::vector<LSHTable<Point>> m_tables;
    size_t m_max_curve_length;
};

#endif  // CURVEDIST_HASH_H
