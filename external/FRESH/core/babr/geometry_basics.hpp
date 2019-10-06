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
#ifndef GEOMETRY_BASICS_HPP_INCLUDED
#define GEOMETRY_BASICS_HPP_INCLUDED

/*
 * Defines the classes point, curve and a function intersection_interval to calculate intersections between points and circles, as well as some helper function and operators.
 */

#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include <cmath>
#include "../types.h"

namespace babr {
 
/*
 * The datatype to use for all calculations that require floating point numbers.
 * double seems to be enough for all testet inputs.
 * For higher precision, long double or boost::multiprecision could be used.
 */
typedef double distance_t;

// bounding box
extern distance_t min_x, min_y, max_x, max_y;

template<typename Point>
Point& operator-=(Point& a, const Point& b);

template<>
inline Point1D& operator-=(Point1D& a, const Point1D& b)
{
    a.x -= b.x;
    return a;
}

template<>
inline Point2D& operator-=(Point2D& a, const Point2D& b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

template<typename Point>
inline Point operator-(const Point& a, const Point& b)
{
    auto tmp = a;
    tmp -= b;
    return tmp;
}

template<typename Point>
Point& operator+=(Point& a, const Point& b);

template<>
inline Point1D& operator+=(Point1D& a, const Point1D& b)
{
    a.x += b.x;
    return a;
}

template<>
inline Point2D& operator+=(Point2D& a, const Point2D& b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

template<typename Point>
inline Point operator+(const Point& a, const Point& b)
{
    auto tmp = a;
    tmp += b;
    return tmp;
}

template<typename Point>
Point& operator/=(Point& a, distance_t d);

template<>
inline Point1D& operator/=(Point1D& a, distance_t d)
{
    a.x /= d;
    return a;
}

template<>
inline Point2D& operator/=(Point2D& a, distance_t d)
{
    a.x /= d;
    a.y /= d;
    return a;
}

inline std::ostream& operator<<(std::ostream& out, const Point1D& p)
{
    out << "(" << p.x << ")";
    return out;
}

inline std::ostream& operator<<(std::ostream& out, const Point2D& p)
{
    out << "(" << p.x << ", " << p.y << ")";
    return out;
}

template<typename T>
inline T sqr(T d)
{
    return d * d;
}

template<typename Point>
distance_t dist(const Point& a, const Point& b);

template<>
inline distance_t dist(const Point1D& a, const Point1D& b)
{
    return std::abs(a.x - b.x);
}

template<>
inline distance_t dist(const Point2D& a, const Point2D& b)
{
    return std::sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
}

template<typename Point>
distance_t dist_sqr(const Point& a, const Point& b);

template<>
inline distance_t dist_sqr(const Point1D& a, const Point1D& b)
{
    return sqr(a.x - b.x);
}

template<>
inline distance_t dist_sqr(const Point2D& a, const Point2D& b)
{
    return sqr(a.x - b.x) + sqr(a.y - b.y);
}

typedef std::pair<distance_t, distance_t> interval; // .first is the startpoint, .second the endpoint (start inclusive, end exclusive)

// This was a constexpr in the original code, but with clang it does not compile
static const interval empty_interval {std::numeric_limits<distance_t>::max(), std::numeric_limits<distance_t>::lowest()};

inline std::ostream& operator<<(std::ostream& out, const interval& i)
{
    out << "[" << i.first << ", " << i.second << "]";
    return out;
}

inline bool is_empty_interval(const interval& i)
{
    return i.first >= i.second;
}

using namespace std;


/*
 * Returns which section of the line segment from line_start to line_end is inside the circle given by circle_center and radius.
 * If the whole line segment lies inside the circle, the result would be [0, 1].
 * If the circle and line segment do not intersect, the result is the empty interval.
 */
template<typename Point>
interval intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end);

template<>
inline interval intersection_interval(Point1D circle_center, distance_t radius, Point1D line_start, Point1D line_end) {
    // move the circle center to (0, 0) to simplify the calculation
    line_start -= circle_center;
    line_end -= circle_center;
    circle_center = {0};

    // The line can be represented as line_start + lambda * v
    const Point1D v = line_end - line_start;

    const distance_t a = sqr(v.x);
    const distance_t b = (line_start.x * v.x);
    const distance_t c = sqr(line_start.x) - sqr(radius);

    const distance_t discriminant = sqr(b / a) - c / a;

    if (discriminant < 0) {
        return empty_interval; // no intersection;
    }

    const distance_t lambda1 = - b / a - sqrt(discriminant);
    const distance_t lambda2 = - b / a + sqrt(discriminant);

    if (lambda2 < 0 || lambda1 > 1) return empty_interval;
    else return {max<distance_t>(lambda1, 0), min<distance_t>(lambda2, 1)};
}

template<>
inline interval intersection_interval(Point2D circle_center, distance_t radius, Point2D line_start, Point2D line_end) {
    // move the circle center to (0, 0) to simplify the calculation
    line_start -= circle_center;
    line_end -= circle_center;
    circle_center = {0, 0};

    // The line can be represented as line_start + lambda * v
    const Point2D v = line_end - line_start;

    // Find Point2Ds p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 = radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2) + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) = radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 * line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, b := 2 * line_start.x * v.x + 2 * line_start.y * v.y, c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * b + c) = 0
    // <=> lambda^2 + (b / a) * lambda + c / a) = 0
    // <=> lambda1/2 = - (b / 2a) +/- sqrt((b / 2a)^2 - c / a)

    const distance_t a = sqr(v.x) + sqr(v.y);
    const distance_t b = (line_start.x * v.x + line_start.y * v.y);
    const distance_t c = sqr(line_start.x) + sqr(line_start.y) - sqr(radius);

    const distance_t discriminant = sqr(b / a) - c / a;

    if (discriminant < 0) {
        return empty_interval; // no intersection;
    }

    const distance_t lambda1 = - b / a - sqrt(discriminant);
    const distance_t lambda2 = - b / a + sqrt(discriminant);

    if (lambda2 < 0 || lambda1 > 1) return empty_interval;
    else return {max<distance_t>(lambda1, 0), min<distance_t>(lambda2, 1)};
}

} // namespace babr

#endif // GEOMETRY_BASICS_HPP_INCLUDED
