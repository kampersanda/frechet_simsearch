#pragma once

#include <FRESH/core/frechet.h>

extern std::vector<double> g_prev_row;
extern std::vector<double> g_cur_row;
#pragma omp threadprivate(g_prev_row, g_cur_row)
std::vector<double> g_prev_row;
std::vector<double> g_cur_row;

template <class Point>
inline bool verify_frechet_continuous(const Curve<Point>& curve_A, const Curve<Point>& curve_B, double range) {
    std::vector<double> dummy;
    return continuous_frechet_distance_predicate(curve_A, curve_B, range, dummy);
}

template <class Point>
inline bool verify_frechet_discrete(const Curve<Point>& curve_A, const Curve<Point>& curve_B, double range) {
    double dist = discrete_frechet_distance(curve_A, curve_B, range, g_prev_row, g_cur_row);
    return 0.0 <= dist && dist <= range;
}

template <class Point>
inline double compute_frechet_discrete(const Curve<Point>& curve_A, const Curve<Point>& curve_B) {
    return discrete_frechet_distance(curve_A.points, curve_B.points, g_prev_row, g_cur_row);
}

template <class Point>
inline double compute_frechet_discrete(const Curve<Point>& curve_A, const Curve<Point>& curve_B, double range) {
    return discrete_frechet_distance(curve_A, curve_B, range, g_prev_row, g_cur_row);
}
