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

#ifndef CURVEDIST_QUERY_H
#define CURVEDIST_QUERY_H

#include "experiment_reporter.h"
#include "frechet.h"
#include "hash.h"
#include "io.h"
#include "query.h"
#include "rand.h"
#include "stats.h"
#include "threading.h"
#include "timer.h"
#include "types.h"

template<typename Point>
inline size_t max_curve_length(const std::vector<Curve<Point>> &data) {
  size_t max_cl = 0;
  for (Curve<Point> c : data) {
    max_cl = std::max(c.points.size(), max_cl);
  }
  return max_cl;
}


/// Positive callback is invoked for each positive match with the index of
/// the query curve as the first argument and the index of the matching dataset curve as the second.
template<typename Point>
void query_continuous_exact(const std::vector<Curve<Point>> &dataset,
                            std::vector<Curve<Point>> &queries,
                            double range,
                            const std::vector<double> &simplification_epsilons,
                            Xorshift1024star &rnd,
                            std::function<void(size_t, size_t)> positive_callback) {
  START_TIMER(algorithm);
  const size_t n = dataset.size();

  for (auto eps : simplification_epsilons) {
    std::cout << " " << eps;
  }
  std::cout << std::endl;

  START_TIMER(all_queries);
  size_t query_idx = 0;
  while (query_idx < queries.size()) {
    const Curve<Point> &query = queries[query_idx];
    START_TIMER(single_query);
    // Process the candidates
#pragma omp parallel for
    for (size_t curve_idx = 0; curve_idx < n; curve_idx++) {
      const Curve<Point> &candidate = dataset[curve_idx];
      if (continuous_frechet_distance_predicate(query, candidate, range, simplification_epsilons)) {
        positive_callback(query_idx, curve_idx);
      }
    }
    query_idx++;
    if (query_idx % 100 == 0) {
      printf("%lu/%lu\n", query_idx, queries.size());
    }
    STOP_TIMER(single_query);
  }
  STOP_TIMER(algorithm);
  STOP_TIMER_V(all_queries);
}

template<typename Point>
void query_continuous_basic(const std::vector<Curve<Point>> &dataset,
                            std::vector<Curve<Point>> &queries,
                            double range,
                            Xorshift1024star &rnd,
                            std::function<void(size_t, size_t)> positive_callback) {
  START_TIMER(algorithm);
  const size_t n = dataset.size();

  START_TIMER(all_queries);
  size_t query_idx = 0;
  while (query_idx < queries.size()) {
    const Curve<Point> &query = queries[query_idx];
    START_TIMER(single_query);
    // Process the candidates
#pragma omp parallel for
    for (size_t curve_idx = 0; curve_idx < n; curve_idx++) {
      const Curve<Point> &candidate = dataset[curve_idx];
      if (babr::get_frechet_distance_upper_bound(query, candidate) <= range){
        positive_callback(query_idx, curve_idx);
      } else if (babr::negfilter(query, candidate, range)) {
        if (query.id == candidate.id) printf("WARNING, ignoring match of %lu with itself (negfilter)", query.id);
        continue;
      } else if (babr::is_frechet_distance_at_most(query, candidate, range)) {
        positive_callback(query_idx, curve_idx);
      } else {
        if (query.id == candidate.id) printf("WARNING, ignoring match of %lu with itself (complete)", query.id);
      }
    }
    query_idx++;
    if (query_idx % 100 == 0) {
      printf("%lu/%lu\n", query_idx, queries.size());
    }
    STOP_TIMER(single_query);
  }
  STOP_TIMER(algorithm);
  STOP_TIMER_V(all_queries);
}

extern std::vector<bool> tl_visited;
#pragma omp threadprivate(tl_visited)
extern std::vector<size_t> tl_collision_counters;
#pragma omp threadprivate(tl_collision_counters)

template<typename Point>
void query_hash(const std::vector<Curve<Point>> &dataset,
                std::vector<Curve<Point>> &queries,
                double range,
                const std::vector<double> &simplification_epsilons,
                size_t k,
                size_t L,
                double resolution_factor,
                Xorshift1024star &rnd,
                std::function<void(size_t, size_t)> positive_callback) {
  START_TIMER(algorithm);
  const size_t n_dataset = dataset.size();
  const size_t max_len = std::max(max_curve_length(dataset), max_curve_length(queries));

#pragma omp parallel
  {
    tl_visited.resize(n_dataset);
  }

  std::cout << "Building table with resolution factor" << resolution_factor << std::endl;
  MultiTable<Point> lsh_table(k, L, resolution_factor * range, n_dataset, max_len, rnd);
  START_TIMER(table_population);
  lsh_table.put_all(dataset);
  STOP_TIMER_V(table_population);

  START_TIMER(all_queries);
#pragma omp parallel for
  for (size_t query_idx = 0; query_idx < queries.size(); query_idx++) {
    const Curve<Point> &query = queries[query_idx];
    START_TIMER(single_query);

    lsh_table.foreach_collision(query.points, tl_visited, [&](size_t collision_idx) {
        if (continuous_frechet_distance_predicate(query, dataset[collision_idx], range,
                                                  simplification_epsilons)) {
          positive_callback(query_idx, collision_idx);
        }
    });

    STOP_TIMER(single_query);
  }
  STOP_TIMER_V(all_queries);
  STOP_TIMER(algorithm);
}


template<typename Point>
void query_hash_no_eval(const std::vector<Curve<Point>> &dataset,
                        std::vector<Curve<Point>> &queries,
                        double range,
                        const std::vector<double> &simplification_epsilons,
                        size_t k,
                        size_t L,
                        double resolution_factor,
                        Xorshift1024star &rnd,
                        std::function<void(size_t, size_t, double)> callback) {
  START_TIMER(algorithm);
  const size_t n_dataset = dataset.size();
  const size_t max_len = std::max(max_curve_length(dataset), max_curve_length(queries));

#pragma omp parallel
  {
    tl_collision_counters.resize(n_dataset);
  }

  std::cout << "Building table with resolution factor" << resolution_factor << std::endl;
  MultiTable<Point> lsh_table(k, L, resolution_factor * range, n_dataset, max_len, rnd);
  START_TIMER(table_population);
  lsh_table.put_all(dataset);
  STOP_TIMER_V(table_population);

  START_TIMER(all_queries);
#pragma omp parallel for
  for (size_t query_idx = 0; query_idx < queries.size(); query_idx++) {
    const Curve<Point> &query = queries[query_idx];
    START_TIMER(single_query);

    lsh_table.count_collisions(query.points, tl_collision_counters);
    for (size_t data_idx=0; data_idx<n_dataset; data_idx++) {
      double score = tl_collision_counters[data_idx] / ((double) L);
      callback(query_idx, data_idx, score);
    }

    /* lsh_table.foreach_collision(query.points, tl_visited, [&](size_t collision_idx) { */
    /*   positive_callback(query_idx, collision_idx); */
    /* }); */

    STOP_TIMER(single_query);
  }
  STOP_TIMER_V(all_queries);
  STOP_TIMER(algorithm);
}


#endif // CURVEDIST_QUERY_H
