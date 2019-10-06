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
#include "frechet.h"
#include <array>
#include "experiment_reporter.h"

// Scratch curves to be used in simplification. They are thread local.
// Reusing them makes the code more cache friendly and reduces allocations
Curve<Point2D> tl_a_simplification_2d;
Curve<Point2D> tl_b_simplification_2d;

Curve<Point1D> tl_a_simplification_1d;
Curve<Point1D> tl_b_simplification_1d;

double c(std::vector<std::vector<double>>& ca, const std::vector<Point2D>& a, const std::vector<Point2D>& b, size_t i,
         size_t j) {
    if (ca[i][j] > -1) {
        return ca[i][j];
    }

    if (i == 0 && j == 0) {
        ca[i][j] = euclidean(a[0], b[0]);
    } else if (i > 0 && j == 0) {
        double rec = c(ca, a, b, i - 1, 0);
        ca[i][j] = std::max(rec, euclidean(a[i], b[0]));
    } else if (i == 0 && j > 0) {
        double rec = c(ca, a, b, 0, j - 1);
        ca[i][j] = std::max(rec, euclidean(a[0], b[j]));
    } else {
        double rec_a = c(ca, a, b, i - 1, j - 1);
        double rec_b = c(ca, a, b, i - 1, j);
        double rec_c = c(ca, a, b, i, j - 1);
        ca[i][j] = std::max(std::min(std::min(rec_a, rec_b), rec_c), euclidean(a[i], b[j]));
    }
    return ca[i][j];
}

double discrete_frechet_distance_r(const std::vector<Point2D>& a, const std::vector<Point2D>& b) {
    size_t n = a.size();
    size_t m = b.size();
    std::vector<std::vector<double>> ca;
    for (size_t i = 0; i < n; i++) {
        std::vector<double> col;
        for (size_t j = 0; j < m; j++) {
            col.push_back(-1);
        }
        ca.push_back(col);
    }
    return c(ca, a, b, n - 1, m - 1);
}

std::vector<size_t> g_tl_stack;
std::vector<bool> g_tl_visited;

std::vector<double> tl_cumulative_dist_a;
std::vector<double> tl_cumulative_dist_b;

int tl_cfdp_current_simplification = 0;
CFDPExitReason tl_cfdp_exit_reason;

size_t tl_cfdp_positives_greedy = 0, tl_cfdp_negatives_filter = 0, tl_cfdp_invocations = 0,
       tl_cfdp_negatives_start_end = 0, tl_cfdp_positives_equal_time = 0, tl_cfdp_negatives_full = 0,
       tl_cfdp_positives_full = 0;

std::array<size_t, MAX_COUNTABLE_SIMPLIFICATIONS>
    tl_cfdp_simplification_positives_equal_time = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    tl_cfdp_simplification_positives_greedy = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    tl_cfdp_simplification_positives_full = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    tl_cfdp_simplification_negatives_filter = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    tl_cfdp_simplification_negatives_full = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

#define _APPEND_CNT_ROW(counter, count) \
    EXPERIMENT_APPEND("continuous-frechet-counters", {{"counter", std::string(counter)}, {"count", count}})

#define _PRINT_CNT_ROW(counter, count) printf("  %s  ::  %lu\n", counter, count)
#define PRINT_CNT_ROW(counter) _PRINT_CNT_ROW(#counter, counter)
#define PRINT_CNT_ARR_ROW(counter, i)                                                                       \
    if (counter[i] != 0) {                                                                                  \
        _PRINT_CNT_ROW((std::to_string(i) + std::string("_") + std::string(#counter)).c_str(), counter[i]); \
    }

#define APPEND_CNT_ROW(counter) _APPEND_CNT_ROW(#counter, counter)
#define APPEND_CNT_ARR_ROW(counter, i)                                                             \
    if (counter[i] != 0) {                                                                         \
        _APPEND_CNT_ROW(std::to_string(i) + std::string("_") + std::string(#counter), counter[i]); \
    }

#define _SUM_TL(variable) variable += tl_cfdp_##variable

void report_counters() {
    size_t positives_greedy = 0, negatives_filter = 0, invocations = 0, negatives_start_end = 0,
           positives_equal_time = 0, negatives_full = 0, positives_full = 0;
    std::array<size_t, MAX_COUNTABLE_SIMPLIFICATIONS>
        simplification_positives_equal_time = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        simplification_positives_greedy = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        simplification_positives_full = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        simplification_negatives_filter = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        simplification_negatives_full = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    size_t checker = 0;

#pragma omp parallel
    {
#pragma omp critical
        {
            _SUM_TL(invocations);
            _SUM_TL(positives_greedy);
            _SUM_TL(negatives_filter);
            _SUM_TL(negatives_start_end);
            _SUM_TL(positives_equal_time);
            _SUM_TL(positives_full);
            _SUM_TL(negatives_full);
            for (size_t i = 0; i < MAX_COUNTABLE_SIMPLIFICATIONS; i++) {
                _SUM_TL(simplification_positives_equal_time[i]);
                _SUM_TL(simplification_positives_greedy[i]);
                _SUM_TL(simplification_positives_full[i]);
                _SUM_TL(simplification_negatives_filter[i]);
                _SUM_TL(simplification_negatives_full[i]);
            }
        }
    }
    APPEND_CNT_ROW(invocations);
    APPEND_CNT_ROW(negatives_start_end);
    checker += negatives_start_end;
    for (size_t i = 0; i < MAX_COUNTABLE_SIMPLIFICATIONS; i++) {
        APPEND_CNT_ARR_ROW(simplification_positives_equal_time, i);
        APPEND_CNT_ARR_ROW(simplification_positives_greedy, i);
        APPEND_CNT_ARR_ROW(simplification_positives_full, i);
        APPEND_CNT_ARR_ROW(simplification_negatives_filter, i);
        APPEND_CNT_ARR_ROW(simplification_negatives_full, i);
        checker += simplification_positives_equal_time[i] + simplification_positives_greedy[i] +
                   simplification_positives_full[i] + simplification_negatives_filter[i] +
                   simplification_negatives_full[i];
    }
    APPEND_CNT_ROW(positives_equal_time);
    APPEND_CNT_ROW(positives_greedy);
    APPEND_CNT_ROW(negatives_filter);
    APPEND_CNT_ROW(negatives_full);
    APPEND_CNT_ROW(positives_full);
    checker += positives_equal_time + positives_greedy + negatives_filter + negatives_full + positives_full;
    if (checker != invocations) {
        // Modified by Shunsuke Kanda (https://github.com/kampersanda)
        printf("WARNING checker != invocations: %zu != %zu\n", checker, invocations);
        // printf("WARNING checker != invocations: %d != %d\n", checker, invocations);
    }

    PRINT_CNT_ROW(invocations);
    PRINT_CNT_ROW(negatives_start_end);
    for (size_t i = 0; i < MAX_COUNTABLE_SIMPLIFICATIONS; i++) {
        PRINT_CNT_ARR_ROW(simplification_positives_equal_time, i);
        PRINT_CNT_ARR_ROW(simplification_positives_greedy, i);
        PRINT_CNT_ARR_ROW(simplification_positives_full, i);
        PRINT_CNT_ARR_ROW(simplification_negatives_filter, i);
        PRINT_CNT_ARR_ROW(simplification_negatives_full, i);
    }
    PRINT_CNT_ROW(positives_equal_time);
    PRINT_CNT_ROW(positives_greedy);
    PRINT_CNT_ROW(negatives_filter);
    PRINT_CNT_ROW(negatives_full);
    PRINT_CNT_ROW(positives_full);
}
