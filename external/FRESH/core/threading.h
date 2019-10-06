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

#ifndef CURVEDIST_THREADING_H
#define CURVEDIST_THREADING_H

#include "types.h"
#include "prelude.h"
#include "rand.h"
#include <omp.h>

struct ThreadState {

  ThreadState(Xorshift1024star rnd)
      : rnd(rnd), scratch_curve(std::vector<Point2D>()),
        frechet_row_1(std::vector<double>()),
        frechet_row_2(std::vector<double>()),
        flags(std::vector<bool>()),
        stack(std::vector<size_t>()) {}

  Xorshift1024star rnd;
  std::vector<Point2D> scratch_curve;
  std::vector<double> frechet_row_1;
  std::vector<double> frechet_row_2;
  std::vector<bool> flags;
  std::vector<size_t> stack;
};

inline std::vector<ThreadState> new_thread_states(Xorshift1024star &rnd) {
  std::vector<ThreadState> thread_states;
  for (int tid = 0; tid < omp_get_max_threads(); tid++) {
    rnd.jump();
    thread_states.emplace_back(ThreadState(rnd));
  }
  return thread_states;
}

#endif // CURVEDIST_THREADING_H
