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

#ifndef CURVEDIST_FRECHER_H
#define CURVEDIST_FRECHER_H

#include "babr/frechet_distance.hpp"
#include "omp.h"
#include "prelude.h"
#include "types.h"

#define INCREASE_COUNTER(counter)                                                                 \
    if (tl_cfdp_current_simplification < 0) {                                                     \
        _Pragma("omp atomic") tl_cfdp_##counter++;                                                \
    } else {                                                                                      \
        _Pragma("omp atomic") tl_cfdp_simplification_##counter[tl_cfdp_current_simplification]++; \
    }

enum CFDPExitReason { GREEDY, NEGATIVE_FILTER, FULL };

extern CFDPExitReason tl_cfdp_exit_reason;
#pragma omp threadprivate(tl_cfdp_exit_reason)

extern int tl_cfdp_current_simplification;
#pragma omp threadprivate(tl_cfdp_current_simplification)
// Thread local counters for frechet distance predicate, for profiling purposes
extern size_t tl_cfdp_positives_greedy, tl_cfdp_negatives_filter, tl_cfdp_invocations, tl_cfdp_negatives_start_end,
    tl_cfdp_positives_equal_time, tl_cfdp_negatives_full, tl_cfdp_positives_full;
#pragma omp threadprivate(tl_cfdp_positives_greedy, tl_cfdp_negatives_filter, tl_cfdp_invocations,           \
                          tl_cfdp_negatives_start_end, tl_cfdp_positives_equal_time, tl_cfdp_negatives_full, \
                          tl_cfdp_positives_full)

#define MAX_COUNTABLE_SIMPLIFICATIONS 16

extern std::array<size_t, MAX_COUNTABLE_SIMPLIFICATIONS> tl_cfdp_simplification_positives_equal_time,
    tl_cfdp_simplification_positives_greedy, tl_cfdp_simplification_positives_full,
    tl_cfdp_simplification_negatives_filter, tl_cfdp_simplification_negatives_full;
#pragma omp threadprivate(tl_cfdp_simplification_positives_equal_time, tl_cfdp_simplification_positives_greedy, \
                          tl_cfdp_simplification_positives_full, tl_cfdp_simplification_negatives_filter,       \
                          tl_cfdp_simplification_negatives_full)

/// Return true if the first or last pair of points are farther than range
template <typename Point>
bool start_end_heuristic(const Curve<Point>& a, const Curve<Point>& b, double range) {
    return ((euclidean(a[0], b[0]) > range) || (euclidean(a[a.size() - 1], b[b.size() - 1]) > range));
}

double discrete_frechet_distance_r(const std::vector<Point2D>& a, const std::vector<Point2D>& b);

template <typename Point>
double discrete_frechet_distance(const std::vector<Point>& a, const std::vector<Point>& b,
                                 std::vector<double>& previous_row, std::vector<double>& current_row) {
    // TODO: see if it's convenient to always use the shortest or
    // longest curve to populate the row

    previous_row.clear();
    current_row.clear();
    previous_row.resize(a.size());
    current_row.resize(a.size());

    // Initialize the first row of the table
    current_row[0] = euclidean_squared(a[0], b[0]);
    for (size_t i = 1; i < a.size(); i++) {
        current_row[i] = std::max(current_row[i - 1], euclidean_squared(a[i], b[0]));
    }

    for (size_t j = 1; j < b.size(); j++) {
        // Swap current and previous row
        current_row.swap(previous_row);

        // Initialize first element of current row
        current_row[0] = std::max(previous_row[0], euclidean_squared(a[0], b[j]));
        for (size_t i = 1; i < a.size(); i++) {
            double min = std::min(current_row[i - 1], std::min(previous_row[i - 1], previous_row[i]));
            current_row[i] = std::max(min, euclidean_squared(a[i], b[j]));
        }
    }

    return std::sqrt(current_row[a.size() - 1]);
}

template <typename Point>
double discrete_frechet_distance(const Curve<Point>& a, const Curve<Point>& b, double bound,
                                 std::vector<double>& previous_row, std::vector<double>& current_row) {
    // TODO: see if it's convenient to always use the shortest or
    // longest curve to populate the row

    // we are using squared Euclidean distances. The additive term is to
    // take into account rounding errors
    bound = (bound * bound) + 0.00001;

    previous_row.clear();
    current_row.clear();
    previous_row.resize(a.size());
    current_row.resize(a.size());

    current_row[0] = euclidean_squared(a[0], b[0]);
    // Check for early termination
    if (current_row[0] > bound || euclidean_squared(a[a.size() - 1], b[b.size() - 1]) > bound) {
        return -1.0;
    }
    // Initialize the first row of the table
    for (size_t i = 1; i < a.size(); i++) {
        current_row[i] = std::max(current_row[i - 1], euclidean_squared(a[i], b[0]));
    }

    for (size_t j = 1; j < b.size(); j++) {
        // Swap current and previous row
        current_row.swap(previous_row);

        // Initialize first element of current row
        current_row[0] = std::max(previous_row[0], euclidean_squared(a[0], b[j]));
        for (size_t i = 1; i < a.size(); i++) {
            double min = std::min(current_row[i - 1], std::min(previous_row[i - 1], previous_row[i]));
            current_row[i] = std::max(min, euclidean_squared(a[i], b[j]));
        }
    }

    return std::sqrt(current_row[a.size() - 1]);
}

inline void idx_to_pair(const size_t k, const size_t wrap, size_t& i, size_t& j) {
    i = k / wrap;
    j = k % wrap;
}

inline void pair_to_idx(const size_t i, const size_t j, const size_t wrap, size_t& k) {
    k = i * wrap + j;
}

template <typename Point>
bool discrete_frechet_distance_predicate(const Curve<Point>& a, const Curve<Point>& b, double bound,
                                         std::vector<size_t>& stack, std::vector<bool>& visited) {
    /* bound = (bound * bound); */
    stack.clear();
    stack.reserve(std::min(a.size(), b.size()));
    visited.clear();
    size_t matrix_size = a.size() * b.size();
    visited.resize(matrix_size);
    std::fill(visited.begin(), visited.end(), false);

    const size_t wrap = b.size();
    size_t i = a.size() - 1;
    size_t j = b.size() - 1;
    if (euclidean(a[i], b[j]) > bound || euclidean(a[0], b[0]) > bound) {
        return false;
    }
    size_t k;
    pair_to_idx(i, j, wrap, k);
    // printf("[%lu] (%lu, %lu) {a: %lu, b: %lu, flags: %lu}\n", k, i, j,
    // a.size(),
    //        b.size(), tl_visited.size());

    visited[k] = true;
    stack.push_back(k);

    while (!stack.empty()) {
        k = stack.back();
        stack.pop_back();
        idx_to_pair(k, wrap, i, j);
        // printf("[%lu] (%lu, %lu) {a: %lu, b: %lu, flags: %lu}\n", k, i, j,
        // a.size(),
        //        b.size(), tl_visited.size());
        if (i == 0 && j == 0) {
            return true;  // We got there!
        }

        if (i > 0 && j > 0) {
            pair_to_idx(i - 1, j - 1, wrap, k);
            if (!visited[k] && euclidean(a[i - 1], b[j - 1]) <= bound) {
                stack.push_back(k);
                visited[k] = true;
            }
        }

        if (i > 0) {
            pair_to_idx(i - 1, j, wrap, k);
            if (!visited[k] && euclidean(a[i - 1], b[j]) <= bound) {
                stack.push_back(k);
                visited[k] = true;
            }
        }

        if (j > 0) {
            pair_to_idx(i, j - 1, wrap, k);
            if (!visited[k] && euclidean(a[i], b[j - 1]) <= bound) {
                stack.push_back(k);
                visited[k] = true;
            }
        }
    }

    // we found no path, fail
    return false;
}

extern std::vector<size_t> g_tl_stack;
extern std::vector<bool> g_tl_visited;
#pragma omp threadprivate(g_tl_stack, g_tl_visited)

template <typename Point>
bool discrete_frechet_distance_predicate(const Curve<Point>& a, const Curve<Point>& b, double bound) {
    return discrete_frechet_distance_predicate(a, b, bound, g_tl_stack, g_tl_visited);
}

/// Report counters to the global experiment
void report_counters();

template <typename Point>
void simplify(const double mu, const Curve<Point>& input, Curve<Point>& output) {
    double mu_mu = mu * mu;
    size_t n = input.size();
    output.clear();
    Point last = input.at(0);
    output.push_back(input[0]);
    for (size_t i = 1; i < n - 1; i++) {
        if (euclidean_squared<Point>(last, input[i]) > mu_mu) {
            output.push_back(input[i]);
            last = input[i];
        }
    }
    output.push_back(input[n - 1]);
}

/// Implementation of the continuous Frechet distance predicate
template <typename Point>
bool continuous_frechet_distance_predicate_impl(const Curve<Point>& a, const Curve<Point>& b, double bound) {
    if (babr::get_frechet_distance_upper_bound(a, b) <= bound) {
        tl_cfdp_exit_reason = CFDPExitReason::GREEDY;
        return true;
    } else if (babr::negfilter(a, b, bound)) {
        tl_cfdp_exit_reason = CFDPExitReason::NEGATIVE_FILTER;
        return false;
    } else if (babr::is_frechet_distance_at_most(a, b, bound)) {
        tl_cfdp_exit_reason = CFDPExitReason::FULL;
        return true;
    } else {
        tl_cfdp_exit_reason = CFDPExitReason::FULL;
        return false;
    }
}

template <typename Point>
bool continuous_frechet_distance_predicate(const Curve<Point>& a, const Curve<Point>& b, double bound);

template <>
inline bool continuous_frechet_distance_predicate(const Curve<Point1D>& a, const Curve<Point1D>& b, double bound) {
    return continuous_frechet_distance_predicate_impl<Point1D>(a, b, bound);
}

template <>
inline bool continuous_frechet_distance_predicate(const Curve<Point2D>& a, const Curve<Point2D>& b, double bound) {
    return continuous_frechet_distance_predicate_impl<Point2D>(a, b, bound);
}

extern std::vector<double> tl_cumulative_dist_a;
extern std::vector<double> tl_cumulative_dist_b;
#pragma omp threadprivate(tl_cumulative_dist_a, tl_cumulative_dist_b)

template <typename Point>
double equal_time_distance(const Curve<Point>& a, const Curve<Point>& b);

template <>
inline double equal_time_distance(const Curve<Point2D>& a, const Curve<Point2D>& b) {
    tl_cumulative_dist_a.clear();
    tl_cumulative_dist_b.clear();
    tl_cumulative_dist_a.resize(a.size());
    tl_cumulative_dist_b.resize(b.size());

    tl_cumulative_dist_a[0] = 0;
    double sum = 0.0;
    for (size_t i = 1; i < a.size(); i++) {
        sum += sum + euclidean(a[i - 1], a[i]);
        tl_cumulative_dist_a[i] = sum;
    }

    tl_cumulative_dist_b[0] = 0;
    sum = 0.0;
    for (size_t i = 1; i < b.size(); i++) {
        sum += sum + euclidean(b[i - 1], b[i]);
        tl_cumulative_dist_b[i] = sum;
    }

    double total_length_a = tl_cumulative_dist_a[a.size() - 1];
    double total_length_b = tl_cumulative_dist_b[b.size() - 1];

    double max_dist = euclidean(a[0], b[0]);

    size_t idx_a = 1;
    size_t idx_b = 1;

    Point2D pt_a, pt_b;

    while ((idx_a < a.size()) && (idx_b < b.size())) {
        // The position in the interval [0, 1] corresponding to the
        // current index for a and b
        double pos_a = tl_cumulative_dist_a[idx_a] / total_length_a;
        double pos_b = tl_cumulative_dist_b[idx_b] / total_length_b;

        if (pos_a < pos_b) {  // the index of a is behind b in terms of equal-time
            pt_a = a[idx_a];
            // we initialize pt_b to a point mid-segment, where the segment is the one
            double time_offset = pos_a - (tl_cumulative_dist_b[idx_b - 1] / total_length_b);
            double diff_x = b[idx_b].x - b[idx_b - 1].x;
            double diff_y = b[idx_b].y - b[idx_b - 1].y;
            pt_b.x = b[idx_b].x + time_offset * diff_x;
            pt_b.y = b[idx_b].y + time_offset * diff_y;
            idx_a++;
        } else {
            pt_b = b[idx_b];
            double time_offset = pos_b - (tl_cumulative_dist_a[idx_a - 1] / total_length_a);
            double diff_x = a[idx_a].x - a[idx_a - 1].x;
            double diff_y = a[idx_a].y - a[idx_a - 1].y;
            pt_a.x = a[idx_a].x + time_offset * diff_x;
            pt_a.y = a[idx_a].y + time_offset * diff_y;
            idx_b++;
        }
        double d = euclidean(pt_a, pt_b);
        if (d > max_dist) {
            max_dist = d;
        }
    }

    double d = euclidean(a[a.size() - 1], b[b.size() - 1]);
    if (d > max_dist) {
        max_dist = d;
    }

    return max_dist;
}

template <>
inline double equal_time_distance(const Curve<Point1D>& a, const Curve<Point1D>& b) {
    tl_cumulative_dist_a.clear();
    tl_cumulative_dist_b.clear();
    tl_cumulative_dist_a.resize(a.size());
    tl_cumulative_dist_b.resize(b.size());

    tl_cumulative_dist_a[0] = 0;
    double sum = 0.0;
    for (size_t i = 1; i < a.size(); i++) {
        sum += sum + euclidean(a[i - 1], a[i]);
        tl_cumulative_dist_a[i] = sum;
    }

    tl_cumulative_dist_b[0] = 0;
    sum = 0.0;
    for (size_t i = 1; i < b.size(); i++) {
        sum += sum + euclidean(b[i - 1], b[i]);
        tl_cumulative_dist_b[i] = sum;
    }

    double total_length_a = tl_cumulative_dist_a[a.size() - 1];
    double total_length_b = tl_cumulative_dist_b[b.size() - 1];

    double max_dist = euclidean(a[0], b[0]);

    size_t idx_a = 1;
    size_t idx_b = 1;

    Point1D pt_a, pt_b;

    while ((idx_a < a.size()) && (idx_b < b.size())) {
        // The position in the interval [0, 1] corresponding to the
        // current index for a and b
        double pos_a = tl_cumulative_dist_a[idx_a] / total_length_a;
        double pos_b = tl_cumulative_dist_b[idx_b] / total_length_b;

        if (pos_a < pos_b) {  // the index of a is behind b in terms of equal-time
            pt_a = a[idx_a];
            // we initialize pt_b to a point mid-segment, where the segment is the one
            double time_offset = pos_a - (tl_cumulative_dist_b[idx_b - 1] / total_length_b);
            double diff_x = b[idx_b].x - b[idx_b - 1].x;
            pt_b.x = b[idx_b].x + time_offset * diff_x;
            idx_a++;
        } else {
            pt_b = b[idx_b];
            double time_offset = pos_b - (tl_cumulative_dist_a[idx_a - 1] / total_length_a);
            double diff_x = a[idx_a].x - a[idx_a - 1].x;
            pt_a.x = a[idx_a].x + time_offset * diff_x;
            idx_b++;
        }
        double d = euclidean(pt_a, pt_b);
        if (d > max_dist) {
            max_dist = d;
        }
    }

    double d = euclidean(a[a.size() - 1], b[b.size() - 1]);
    if (d > max_dist) {
        max_dist = d;
    }

    return max_dist;
}

enum FuzzyResult { YES, YES_ET, NO, MAYBE };

extern Curve<Point2D> tl_a_simplification_2d;
extern Curve<Point2D> tl_b_simplification_2d;
extern Curve<Point1D> tl_a_simplification_1d;
extern Curve<Point1D> tl_b_simplification_1d;
#pragma omp threadprivate(tl_a_simplification_2d, tl_b_simplification_2d, tl_a_simplification_1d, \
                          tl_b_simplification_1d)

template <typename Point>
FuzzyResult fuzzy_decide(const Curve<Point>& a, const Curve<Point>& b, const double epsilon, const double delta,
                         Curve<Point>& simplification_a, Curve<Point>& simplification_b) {
    double mu = epsilon / 4 * delta;
    double delta_prime = delta + 2 * mu;
    simplify<Point>(mu, a, simplification_a);
    simplify<Point>(mu, b, simplification_b);
    if (!continuous_frechet_distance_predicate(simplification_a, simplification_b, delta_prime)) {
        return FuzzyResult::NO;
    }
    double delta_positive = delta / (1 + epsilon);
    mu = epsilon / 4 * delta_positive;
    delta_prime = delta_positive + 2 * mu;
    simplify<Point>(mu, a, simplification_a);
    simplify<Point>(mu, b, simplification_b);
    if (equal_time_distance(simplification_a, simplification_b) <= delta_prime) {
        return FuzzyResult::YES_ET;
    } else if (continuous_frechet_distance_predicate(simplification_a, simplification_b, delta_prime)) {
        return FuzzyResult::YES;
    } else {
        return FuzzyResult::MAYBE;
    }
}

template <typename Point>
FuzzyResult fuzzy_decide(const Curve<Point>& a, const Curve<Point>& b, const double epsilon, const double delta);

template <>
inline FuzzyResult fuzzy_decide(const Curve<Point1D>& a, const Curve<Point1D>& b, const double epsilon,
                                const double delta) {
    return fuzzy_decide(a, b, epsilon, delta, tl_a_simplification_1d, tl_b_simplification_1d);
}

template <>
inline FuzzyResult fuzzy_decide(const Curve<Point2D>& a, const Curve<Point2D>& b, const double epsilon,
                                const double delta) {
    return fuzzy_decide(a, b, epsilon, delta, tl_a_simplification_2d, tl_b_simplification_2d);
}

template <typename Point>
bool continuous_frechet_distance_predicate(const Curve<Point>& a, const Curve<Point>& b, double range,
                                           const std::vector<double>& simplification_epsilons) {
    const size_t n_simplifications = simplification_epsilons.size();
#pragma omp atomic
    tl_cfdp_invocations++;
    if (start_end_heuristic(a, b, range)) {
#pragma omp atomic
        tl_cfdp_negatives_start_end++;
        return false;  // This is a negative match
    }

    // Check the simplifications
    for (size_t i = 0; i < n_simplifications; i++) {
        tl_cfdp_current_simplification = i;
        double epsilon = simplification_epsilons[i];
        switch (fuzzy_decide(a, b, epsilon, range)) {
            case FuzzyResult::YES_ET:
                INCREASE_COUNTER(positives_equal_time);
                /* tl_cfdp_simplification_positives_equal_time.at(i)++; */
                return true;
            case FuzzyResult::YES:
                /* tl_cfdp_simplification_positives.at(i)++; */
                switch (tl_cfdp_exit_reason) {
                    case CFDPExitReason::GREEDY:
                        INCREASE_COUNTER(positives_greedy);
                        break;
                    case CFDPExitReason::FULL:
                        INCREASE_COUNTER(positives_full);
                        break;
                    default:
                        throw std::logic_error("should not be here");
                }
                return true;
            case FuzzyResult::NO:
                /* tl_cfdp_simplification_negatives.at(i)++; */
                switch (tl_cfdp_exit_reason) {
                    case CFDPExitReason::NEGATIVE_FILTER:
                        INCREASE_COUNTER(negatives_filter);
                        break;
                    case CFDPExitReason::FULL:
                        INCREASE_COUNTER(negatives_full);
                        break;
                    default:
                        throw std::logic_error("should not be here");
                }
                return false;
            case FuzzyResult::MAYBE:
                // Go on and check next simplification
                break;
        }
    }
    tl_cfdp_current_simplification = -1;

    // Check equal time distance
    if (equal_time_distance(a, b) <= range) {
        INCREASE_COUNTER(positives_equal_time);
        /* tl_cfdp_equal_time_positives++; */
        return true;
    }

    // Check the true distance
    if (continuous_frechet_distance_predicate(a, b, range)) {
        /* tl_cfdp_full_distance_positives++; */
        switch (tl_cfdp_exit_reason) {
            case CFDPExitReason::GREEDY:
                INCREASE_COUNTER(positives_greedy);
                break;
            case CFDPExitReason::FULL:
                INCREASE_COUNTER(positives_full);
                break;
            default:
                throw std::logic_error("should not be here");
        }
        return true;
    } else {
        /* tl_cfdp_full_distance_negatives++; */
        switch (tl_cfdp_exit_reason) {
            case CFDPExitReason::NEGATIVE_FILTER:
                INCREASE_COUNTER(negatives_filter);
                break;
            case CFDPExitReason::FULL:
                INCREASE_COUNTER(negatives_full);
                break;
            default:
                throw std::logic_error("should not be here");
        }
        return false;
    }
}

#endif  // CURVEDIST_FRECHER_H
