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
#include "types.h"
#include "threading.h"

QueryResult::QueryResult()
  : query(nullptr), ignored(0),
    matching_curves(std::vector<matching_curves_t>(omp_get_max_threads())) {}

QueryResult::QueryResult(Curve<Point2D> const * query)
  : query(query), ignored(0), matching_curves(std::vector<matching_curves_t>(omp_get_max_threads())) {}

void QueryResult::append(Curve<Point2D> const *c, size_t num_collisions) {
  size_t tid = omp_get_thread_num();
  matching_curves[tid].push_back(std::make_pair(c, num_collisions));
}

std::ostream& operator<<(std::ostream& os, const Point2D &p) {
  os << "(" << p.x << ", " << p.y << ")";
  return os;
}


