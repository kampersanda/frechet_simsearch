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
#include "densify.h"
#include "frechet.h"

void densify(Curve<Point2D> &curve, double max_gap, std::vector<Point2D> &scratch_curve) {
  // first, count the number of points that have to be added, so to be able to
  // directly allocate the space we need without expensive guessing
  size_t points_to_add = 0;
  const size_t n = curve.points.size();
  for (size_t i = 0; i < n - 1; i++) {
    const double d = euclidean(curve.points[i], curve.points[i + 1]);
    if (d > max_gap) {
      // count the points we have to add to fill this gap
      points_to_add += std::floor(d / max_gap);
    }
  }
  scratch_curve.resize(n); // vector capacity is not reduced
  std::copy(curve.points.begin(), curve.points.end(), scratch_curve.begin());
  const size_t new_size = n + points_to_add;
  curve.points.clear();
  curve.points.reserve(new_size);

  for (size_t orig_idx = 0; orig_idx < n - 1; orig_idx++) {
    curve.points.push_back(scratch_curve[orig_idx]);
    const double d =
        euclidean(scratch_curve[orig_idx], scratch_curve[orig_idx + 1]);
    if (d > max_gap) {
      // Add points until the gap is filled
      const size_t n_pts = std::floor(d / max_gap);
      const double base_x = scratch_curve[orig_idx].x;
      const double base_y = scratch_curve[orig_idx].y;
      const double diff_x =
          scratch_curve[orig_idx + 1].x - scratch_curve[orig_idx].x;
      const double diff_y =
          scratch_curve[orig_idx + 1].y - scratch_curve[orig_idx].y;
      for (size_t j = 0; j < n_pts; j++) {
        double x = base_x + (j * max_gap) * diff_x / d;
        double y = base_y + (j * max_gap) * diff_y / d;
        curve.points.emplace_back(x, y);
      }
    }
  }
}
