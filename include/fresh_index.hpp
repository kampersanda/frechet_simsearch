#pragma once

#include <FRESH/core/query.h>

#include "basic.hpp"

template <class Point>
class FreshIndex {
  public:
    using point_type = Point;
    using table_type = MultiTable<Point>;

  public:
    FreshIndex() = default;

    FreshIndex(const std::vector<Curve<Point>>& curves, const Configure& cfg) : m_size(curves.size()) {
        Splitmix64 seeder(cfg.seed);
        Xorshift1024star rnd(seeder.next());
        m_table = std::make_unique<table_type>(cfg.concatenations, cfg.repeatations, cfg.resolution,  //
                                               curves.size(), max_curve_length(curves), rnd);
        m_table->put_all(curves);
    }

    ~FreshIndex() = default;

    void count_collisions(const Curve<Point>& query_curve, std::vector<size_t>& counters) {
        m_table->count_collisions(query_curve.points, counters);
    }

    size_t size() const {
        return m_size;
    }
    size_t get_memory_usage() const {
        return 0;
    }

  private:
    size_t m_size = 0;
    std::unique_ptr<table_type> m_table;
};