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

#ifndef CURVEDIST_TYPES_H
#define CURVEDIST_TYPES_H

#include "msgpack.hpp"
#include "prelude.h"

typedef double distance_t;

struct Point1D {

  static size_t dimensions() { return 1; }

  Point1D() : x(0) {}

  Point1D(double x) : x(x) {}

  double x;

  MSGPACK_DEFINE(x);
};


struct Point2D {

  static size_t dimensions() { return 2; }

  Point2D() : x(0), y(0) {}

  Point2D(double x, double y) : x(x), y(y) {}

  double x;
  double y;

  MSGPACK_DEFINE(x, y);
};

std::ostream& operator<<(std::ostream& os, const Point2D &p);

template<typename Point>
double euclidean(const Point & a, const Point& b);

template<typename Point>
double euclidean_squared(const Point & a, const Point& b);

template<>
inline double euclidean(const Point2D &a, const Point2D &b) {
  double diff_x = a.x - b.x;
  double diff_y = a.y - b.y;
  return std::sqrt(diff_x * diff_x + diff_y * diff_y);
}

template<>
inline double euclidean_squared(const Point2D &a, const Point2D &b) {
  double diff_x = a.x - b.x;
  double diff_y = a.y - b.y;
  return diff_x * diff_x + diff_y * diff_y;
}

template<>
inline double euclidean(const Point1D &a, const Point1D &b) {
  return std::abs(a.x - b.x);
}

template<>
inline double euclidean_squared(const Point1D &a, const Point1D &b) {
  double diff_x = a.x - b.x;
  return diff_x * diff_x;
}


template <typename Point>
struct Curve {

  Curve() : id(0u), points(std::vector<Point>()) {}
  Curve(size_t id)
      : id(id), points(std::vector<Point>()) {}
  Curve(size_t id, const std::vector<Point> &pts)
      : id(id), points(pts) {}

  // I would have liked to have this field declared as const, however this makes
  // the object non-assignable, which doesn't play well with STL containers. For
  // instance one could not sort a vector of Curve structs using std::sort if
  // id were const.
  size_t id;
  std::vector<Point> points;

  MSGPACK_DEFINE_MAP(id, points);

  // to be called right after deserialization from MSGPACK
  void fix_prefix_lengths() {
    prefix_length.resize(points.size());
    prefix_length[0] = 0;
    for (int i = 1; i < points.size(); ++i) {
      prefix_length[i] = prefix_length[i - 1] + euclidean(points[i - 1], points[i]);
    }
  }

  inline size_t size() const { return points.size(); }

  inline Point operator[](size_t i) const
  {
    return points[i];
  }

  inline Point at(size_t i) const
  {
    return points.at(i);
  }

  inline void clear() {
    points.clear();
    prefix_length.clear();
  }

  inline distance_t curve_length(size_t i, size_t j) const
  {
    return prefix_length[j] - prefix_length[i];
  }

  inline Point front() const { return points.front(); }

  inline Point back() const { return points.back(); }

  inline void pop_back() {
    points.pop_back();
    prefix_length.pop_back();
  }

  inline void push_back(Point p) {
    if (prefix_length.size()) {
      prefix_length.push_back(prefix_length.back() + euclidean(points.back(), p));
    } else {
      prefix_length.push_back(0);
    }
    points.push_back(p);
  }


private:
    std::vector<double> prefix_length;

    friend typename std::vector<Point>::const_iterator begin(const Curve<Point>& c);
    friend typename std::vector<Point>::const_iterator end(const Curve<Point>& c);
};

/* template<typename Point> */
/* typename std::vector<Point>::const_iterator begin(const Curve<Point>& c) { */
/*     return c.points.begin(); */
/* } */

/* template<typename Point> */
/* typename std::vector<Point>::const_iterator end(const Curve<Point>& c) { */
/*     return c.points.end(); */
/* } */

inline std::vector<Point1D>::const_iterator begin(const Curve<Point1D>& c) {
    return c.points.begin();
}
inline std::vector<Point1D>::const_iterator end(const Curve<Point1D>& c) {
    return c.points.end();
}

inline std::vector<Point2D>::const_iterator begin(const Curve<Point2D>& c) {
    return c.points.begin();
}
inline std::vector<Point2D>::const_iterator end(const Curve<Point2D>& c) {
    return c.points.end();
}


/* template<typename Point> */
/* inline std::ostream& operator<<(std::ostream& out, const Curve<Point>& c) { */
/*     out << "["; */
/*     for (auto p: c) out << p << ", "; */
/*     out << "]"; */
/*     return out; */
/* } */


struct QueryResult {

  typedef std::vector<std::pair<Curve<Point2D> const *, size_t>> matching_curves_t;

  QueryResult();
  QueryResult(Curve<Point2D> const * query);
  
  void append(Curve<Point2D> const *c, size_t num_collisions);

  void set_time(double t) { time_ms = t; }
  
  Curve<Point2D> const * query;
  size_t ignored;
  double time_ms;
  std::vector<matching_curves_t> matching_curves;

  template<typename Stream>
  void pack(Stream& out) {
    msgpack::packer<Stream> pk(out);
    pk.pack_map(4);
    pk.pack(std::string("query curve"));
    pk.pack_unsigned_long_long(query->id);

    pk.pack(std::string("query ignored curves"));
    if (ignored) pk.pack_true();
    else pk.pack_false();

    pk.pack(std::string("time"));
    pk.pack_double(time_ms);
  
    pk.pack(std::string("count_collisions curves"));
    size_t num_matches = 0;
    for (auto &collection : matching_curves) {
      num_matches += collection.size();
    }
    pk.pack_array(num_matches);
    for (auto &collection : matching_curves) {
      for (auto &pair : collection) {
        pk.pack_array(2);
        pk.pack_unsigned_long_long(pair.first->id);
        pk.pack_unsigned_long_long(pair.second);
      }
    }
  }
};

#endif // CURVEDIST_TYPES_H

