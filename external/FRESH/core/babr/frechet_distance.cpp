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
#include "frechet_distance.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include <cassert>
#include <cstddef>

namespace babr {

using namespace babr;
using std::vector;
using std::numeric_limits;
using std::max;
using std::min;

template<>
distance_t get_frechet_distance(const Curve<Point1D>& a, const Curve<Point1D>& b)
{
    distance_t min_coordinate = a[0].x, max_coordinate = a[0].x;
    for (const Curve<Point1D>& c: {a, b}) {
        for (const Point1D& p: c) {
            min_coordinate = min(min_coordinate, p.x);
            max_coordinate = max(max_coordinate, p.x);
        }
    }

    distance_t min_d = 0, max_d = (max_coordinate - min_coordinate) * 3 + 1;
    int cnt = 0;
    while (min_d + epsilon < max_d) {
        ++cnt;
        distance_t m = (min_d + max_d) / 2;
        if (is_frechet_distance_at_most(a, b, m)) max_d = m;
        else min_d = m;
    }

    return min_d;
}

template<>
distance_t get_frechet_distance(const Curve<Point2D>& a, const Curve<Point2D>& b)
{
    distance_t min_coordinate = a[0].x, max_coordinate = a[0].x;
    for (const Curve<Point2D>& c: {a, b}) {
        for (const Point2D& p: c) {
            min_coordinate = min(min_coordinate, min(p.x, p.y));
            max_coordinate = max(max_coordinate, max(p.x, p.y));
        }
    }

    distance_t min_d = 0, max_d = (max_coordinate - min_coordinate) * 3 + 1;
    int cnt = 0;
    while (min_d + epsilon < max_d) {
        ++cnt;
        distance_t m = (min_d + max_d) / 2;
        if (is_frechet_distance_at_most(a, b, m)) max_d = m;
        else min_d = m;
    }

    return min_d;
}

/* distance_t get_frechet_distance_upper_bound(const curve& a, const curve& b) */
/* { */
/*     distance_t distance = dist(a.back(), b.back()); */
/*     size_t pos_a = 0, pos_b = 0; */

/*     while (pos_a + pos_b < a.size() + b.size() - 2) { */
/*         distance = max(distance, dist(a[pos_a], b[pos_b])); */
/*         if (pos_a == a.size() - 1) ++pos_b; */
/*         else if (pos_b == b.size() - 1) ++pos_a; */
/*         else { */
/*             distance_t dist_a = dist_sqr(a[pos_a + 1], b[pos_b]); */
/*             distance_t dist_b = dist_sqr(a[pos_a], b[pos_b + 1]); */
/*             distance_t dist_both = dist_sqr(a[pos_a + 1], b[pos_b + 1]); */
/*             if (dist_a < dist_b && dist_a < dist_both) { */
/*                 ++pos_a; */
/*             } else if (dist_b < dist_both) { */
/*                 ++pos_b; */
/*             } else { */
/*                 ++pos_a; */
/*                 ++pos_b; */
/*             } */
/*         } */
/*     } */

/*     return distance; */
/* } */

/* size_t nextclosepoint(const curve& c, size_t i, point p, distance_t d) */
/* { */
/* 	size_t delta = 1; */
/* 	size_t k = i; */
/* 	while (true) { */
/* 		if (k == c.size()-1) { */
/* 			if (dist(c[k],p) <= d) { return k; } */
/* 			else { return c.size(); } */
/* 		} */
/* 		else { */
/* 			delta = std::min(delta, c.size() - 1 - k); */
/* 			if (dist(p, c[k]) - c.curve_length(k,k+delta) > d) { */
/* 				k += delta; */
/* 				delta *= 2; */
/* 			} */
/* 			else if (delta > 1) {delta /= 2; } */
/* 			else { return k; } */
/* 		} */
/* 	} */
/* } */

/* bool negfilter(const curve& c1, const curve& c2, distance_t d) */
/* { */
/* 	for (size_t delta = std::max(c1.size(),c2.size())-1; delta >= 1; delta /= 2) { */
/* 		size_t i = 0; */
/* 		for (size_t j = 0; j < c2.size(); j += delta) { */
/* 			i = nextclosepoint(c1, i, c2[j], d); */
/* 			if (i >= c1.size()) { */
/* 				return true; */
/* 			} */
/* 		} */
/* 		size_t j = 0; */
/* 		for (size_t i = 0; i < c1.size(); i += delta) { */
/* 			j = nextclosepoint(c2, j, c1[i], d); */
/* 			if (j >= c2.size()) { */
/* 				return true; */
/* 			} */
/* 		} */
/* 	} */
/* 	return false; */
/* } */

} // namespace babr







