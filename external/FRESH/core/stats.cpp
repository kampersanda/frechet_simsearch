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
#include "stats.h"
#include "frechet.h"
#include "threading.h"

DatasetStats compute_stats(const std::vector<Curve<Point2D>> &dataset,
                           size_t num_sampled_pairs, Xorshift1024star &rnd) {
  DatasetStats stats;
  size_t n = dataset.size();
  size_t total_pairs = n * (n - 1) / 2;
  num_sampled_pairs = std::min(total_pairs, num_sampled_pairs);
  printf("Computing stats on %lu pairs\n", num_sampled_pairs);

  std::vector<double> distances(num_sampled_pairs, -1.0);

  std::vector<ThreadState> t_state = new_thread_states(rnd);

#pragma omp parallel for
  for (size_t h = 0; h < num_sampled_pairs; h++) {
    size_t tid = omp_get_thread_num();
    ThreadState &ts = t_state[tid];
    size_t i = std::floor(ts.rnd.next_double(n));
    size_t j = std::floor(ts.rnd.next_double(n));
    const std::vector<Point2D> &a = dataset[i].points;
    const std::vector<Point2D> &b = dataset[j].points;
    double d = discrete_frechet_distance(a, b, ts.frechet_row_1, ts.frechet_row_2);
    distances[h] = d;
  }

  std::sort(distances.begin(), distances.end());
  stats.max_distance = 0.0;
  stats.min_distance = std::numeric_limits<double>::infinity();
  stats.mean_distance = 0.0;
  stats.median_distance = distances[num_sampled_pairs / 2];
  for (double d : distances) {
    stats.max_distance = std::max(d, stats.max_distance);
    stats.min_distance = std::min(d, stats.min_distance);
    stats.mean_distance += d;
  }
  stats.mean_distance /= num_sampled_pairs;

  return stats;
}
