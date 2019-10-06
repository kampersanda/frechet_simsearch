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
#ifndef FRECHET_DISTANCE_HPP_INCLUDED
#define FRECHET_DISTANCE_HPP_INCLUDED

#include "geometry_basics.hpp"
#include "frechet_distance2.hpp"

namespace babr {
  using namespace babr;

// for get_frechet_distance
constexpr distance_t epsilon = 1e-5l;

/*
 * Checks whether the frechet distance is at most d, with the maximal accuracy distance_t allows.
 * O(a.size() * b.size())
 */
template<typename Point>
bool is_frechet_distance_at_most(const Curve<Point>& a, const Curve<Point>& b, distance_t d) {
    assert(a.size());
    assert(b.size());
    if (dist(a.front(), b.front()) > d || dist(a.back(), b.back()) > d) return false;
    if (a.size() == 1 && b.size() == 1) return true;
    else if (a.size() == 1) return get_dist_to_point_sqr(b, a[0]) <= sqr(d);
    else if (b.size() == 1) return get_dist_to_point_sqr(a, b[0]) <= sqr(d);

    vector<interval> ra, rb, ra_out, rb_out;
    ra.push_back({0, get_last_reachable_point_from_start(a, b, d)});
    rb.push_back({0, get_last_reachable_point_from_start(b, a, d)});

    get_reachable_intervals(0, a.size() - 1, 0, b.size() - 1, a, b, d, ra, rb, ra_out, rb_out);

    return ra_out.size() && (ra_out.back().second >= b.size() - static_cast<distance_t>(1.5));
}
/*
 * Returns the frechet distance between a and b, accurate to +/- epsilon.
 * O(log(possible results) * a.size() * b.size())
 */
template<typename Point>
distance_t get_frechet_distance(const Curve<Point>& a, const Curve<Point>& b);

template<>
distance_t get_frechet_distance(const Curve<Point2D>& a, const Curve<Point2D>& b);
template<>
distance_t get_frechet_distance(const Curve<Point1D>& a, const Curve<Point1D>& b);

/*
 * Calculates an upper bound for the frechet distance of a and b by guessing a matching between a and b
 * O(a.size() + b.size())
 */
template<typename Point>
distance_t get_frechet_distance_upper_bound(const Curve<Point>& a, const Curve<Point>& b) {
    distance_t distance = dist(a.back(), b.back());
    size_t pos_a = 0, pos_b = 0;

    while (pos_a + pos_b < a.size() + b.size() - 2) {
        distance = max(distance, dist(a[pos_a], b[pos_b]));
        if (pos_a == a.size() - 1) ++pos_b;
        else if (pos_b == b.size() - 1) ++pos_a;
        else {
            distance_t dist_a = dist_sqr(a[pos_a + 1], b[pos_b]);
            distance_t dist_b = dist_sqr(a[pos_a], b[pos_b + 1]);
            distance_t dist_both = dist_sqr(a[pos_a + 1], b[pos_b + 1]);
            if (dist_a < dist_b && dist_a < dist_both) {
                ++pos_a;
            } else if (dist_b < dist_both) {
                ++pos_b;
            } else {
                ++pos_a;
                ++pos_b;
            }
        }
    }

    return distance;
}
/*
 * Returns (a discrete approximation of) the first point c[j] on c, with j >= i, that is within distance d of point p.
 */
template<typename Point>
size_t nextclosepoint(const Curve<Point>& c, size_t i, Point p, distance_t d) {
	size_t delta = 1;
	size_t k = i;
	while (true) {
		if (k == c.size()-1) {
			if (dist(c[k],p) <= d) { return k; }
			else { return c.size(); }
		}
		else {
			delta = std::min(delta, c.size() - 1 - k);
			if (dist(p, c[k]) - c.curve_length(k,k+delta) > d) {
				k += delta;
				delta *= 2;
			}
			else if (delta > 1) {delta /= 2; }
			else { return k; }
		}
	}
}

/*
 * Tries to show that the Frechet distance of a and b is more than d. Returns true if a proof is found.
 */
template<typename Point>
bool negfilter(const Curve<Point>& c1, const Curve<Point>& c2, distance_t d) {
	for (size_t delta = std::max(c1.size(),c2.size())-1; delta >= 1; delta /= 2) {
		size_t i = 0;
		for (size_t j = 0; j < c2.size(); j += delta) {
			i = nextclosepoint(c1, i, c2[j], d);
			if (i >= c1.size()) {
				return true;
			}
		}
		size_t j = 0;
		for (size_t i = 0; i < c1.size(); i += delta) {
			j = nextclosepoint(c2, j, c1[i], d);
			if (j >= c2.size()) {
				return true;
			}
		}
	}
	return false;
}

}

#endif // FRECHET_DISTANCE_HPP_INCLUDED
