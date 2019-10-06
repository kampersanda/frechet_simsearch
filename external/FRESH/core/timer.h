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

#ifndef CURVEDIST_TIMER_H
#define CURVEDIST_TIMER_H

#include "prelude.h"
#include <boost/timer/timer.hpp>

double get_milliseconds(const boost::timer::cpu_timer &timer);

class TimeRegistry {

private:
  typedef std::list<boost::timer::cpu_times> elapsed_list_t;
  typedef std::unordered_map<std::string, elapsed_list_t> timers_map_t;

private:
  TimeRegistry() : m_timers(timers_map_t()) {}

public:
  TimeRegistry(TimeRegistry const &) = delete;
  void operator=(TimeRegistry const &) = delete;

  static TimeRegistry &getInstance() {
    // Guaranteed to be destroyed.
    // Instantiated on first use.
    static TimeRegistry instance;
    return instance;
  }

  void append(const std::string &name, boost::timer::cpu_times elapsed) {
    if (m_timers.count(name)) {
      m_timers[name].push_back(elapsed);
    } else {
      elapsed_list_t l;
      l.push_back(elapsed);
      m_timers.insert(std::make_pair(name, l));
    }
  }

  void print_all() const {
    for (const auto &entry : m_timers) {
      printf("=== %s ============\n", entry.first.c_str());
      size_t i = 0;
      for (const auto &t : entry.second) {
        double e = t.wall / 1000000.0;
        printf("  %lu : %f\n", i, e);
        i++;
      }
    }
  }

  std::list<std::pair<std::string, std::list<double>>> get_walls() const {
    std::list<std::pair<std::string, std::list<double>>> walls;
    for (const auto &entry : m_timers) {
      std::list<double> ws;
      std::string name = entry.first;
      size_t i = 0;
      for (const auto &t : entry.second) {
        double e = t.wall / 1000000.0;
        ws.push_back(e);
        i++;
      }
      walls.push_back(std::make_pair(name, ws));
    }
    return walls;
  }

  double get_wall_sum(std::string name) const {
    double sum = 0;
    for (const auto &t : m_timers.at(name)) {
        double wall_time = t.wall / 1000000.0;
        sum += wall_time;
    }
    return sum;
  }

private:
  timers_map_t m_timers;
};

#define START_TIMER(___name) boost::timer::cpu_timer ___name;

#define STOP_TIMER(___name)                                             \
  do {                                                                  \
    ___name.stop();                                                     \
    TimeRegistry::getInstance().append(#___name, ___name.elapsed());    \
  } while (false);

#define GET_TIMER_MS(___name)                      \
  ___name.elapsed().wall / 1000000.0

#define STOP_TIMER_V(___name)                                           \
  do {                                                                  \
    ___name.stop();                                                     \
    TimeRegistry::getInstance().append(#___name, ___name.elapsed());    \
    printf("[%s] %s", #___name, ___name.format(4).c_str());             \
  } while (false);

#define GET_WALL_SUM(name_str)    \
    TimeRegistry::getInstance().get_wall_sum(name_str)


#define EXPERIMENT_RECORD_PROFILE()                                            \
  for (const auto &wall : TimeRegistry::getInstance().get_walls()) {           \
    std::string timer = wall.first;                                            \
    size_t i = 0;                                                              \
    for (const auto &time : wall.second) {                                     \
      EXPERIMENT_APPEND("profile",                                             \
                        {{"timer", timer}, {"measure", i}, {"time", time}});   \
      i++;                                                                     \
    }                                                                          \
  }

#endif // CURVEDIST_TIMER_H
