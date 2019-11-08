#include <cmath>
#include <ctime>
#include <random>

#include <FRESH/core/experiment_reporter.h>
#include <FRESH/core/scores.h>
#include <FRESH/core/timer.h>

#include <cmdline.h>

#include <groundtruth.hpp>
#include <multi_trie_index.hpp>
#include <vcode_array.hpp>

// #define ENABLE_VERIFICATION

using namespace frechet_simsearch;

const double ABORT_THRESHOLD_MS = 100.0;

std::string make_result_path(std::string result_dir, std::string base_path, std::string query_path,
                             double frechet_range, std::string hamming_ranges, Configure cfg) {
    auto result_path = tfm::format("%s/mtrie_perf-%s-%s-%g-%s-%s.mpk",  //
                                   result_dir, normalize_filepath(base_path), normalize_filepath(query_path),  //
                                   frechet_range, normalize_filepath(hamming_ranges), cfg.to_str());
#ifdef EXPERIMENT_COMPRESS
    result_path += ".bz2";
#endif
    return result_path;
}

template <class Point, uint32_t Length>
int main_template(const cmdline::parser& p) {
    tfm::printfln("num threads: %d", omp_get_max_threads());

    const auto base_path = p.get<std::string>("base_path");
    const auto query_path = p.get<std::string>("query_path");
    const auto result_dir = p.get<std::string>("result_dir");
    const auto frechet_range = p.get<double>("frechet_range");
    const auto hamming_ranges = p.get<std::string>("hamming_ranges");
    const auto repeatations = p.get<uint32_t>("repeatations");
    const auto buckets = p.get<uint32_t>("buckets");
    const auto reduction_threshold = p.get<uint32_t>("reduction_threshold");
    const auto resolution_factor = p.get<double>("resolution_factor");
    const auto seed = p.get<uint64_t>("seed");

    const double resolution = resolution_factor * frechet_range * Point::dimensions();

    tfm::printfln("frechet_range: %g", frechet_range);
    tfm::printfln("resolution: %g (factor: %g)", resolution, resolution_factor);

    Configure cfg;
    {
        cfg.repeatations = repeatations;
        cfg.buckets = buckets;
        cfg.reduction_threshold = reduction_threshold;
        cfg.resolution = resolution;
        cfg.resolution_factor = resolution_factor;
        cfg.seed = seed;
    }
    cfg.show(std::cout);

    const std::string result_path = make_result_path(result_dir, base_path, query_path,  //
                                                     frechet_range, hamming_ranges, cfg);
    if (exists_path(result_path)) {
        warnfln("already exists: %s", result_path);
        return 1;
    }

    /**
     *  Load database
     */
    const std::vector<Curve<Point>> base_curves = load_dataset_msgpack<Point>(base_path);
    const std::vector<Curve<Point>> query_curves = load_dataset_msgpack<Point>(query_path);

    const size_t base_size = base_curves.size();
    const size_t query_size = query_curves.size();

    tfm::printfln("number of base curves : %d", base_size);
    tfm::printfln("number of query curves: %d", query_size);

    EXPERIMENT_TAG("algorithm", std::string("MTRIE"));
    EXPERIMENT_TAG("base_path", base_path);
    EXPERIMENT_TAG("query_path", query_path);
    EXPERIMENT_TAG("base_size", base_size);
    EXPERIMENT_TAG("query_size", query_size);
    EXPERIMENT_TAG("frechet_range", frechet_range);
    EXPERIMENT_TAG("hamming_ranges", hamming_ranges);
    EXPERIMENT_TAG("repeatations", size_t(repeatations));
    EXPERIMENT_TAG("buckets", size_t(buckets));
    EXPERIMENT_TAG("reduction_threshold", size_t(reduction_threshold));
    EXPERIMENT_TAG("resolution_factor", resolution_factor);
    EXPERIMENT_TAG("resolution", resolution);
    EXPERIMENT_TAG("seed", size_t(seed));
    EXPERIMENT_TAG("dimensions", Point::dimensions());

    /**
     *  Build index
     */
    START_TIMER(preprocessing);
    MultiTrieIndex<Point, Length> index(base_curves, cfg);
    STOP_TIMER_V(preprocessing);

    const double constr_sec = GET_TIMER_MS(preprocessing) / 1000.0;
    const size_t memory_usage = index.get_memory_usage();
    const size_t index_memory_usage = index.get_index_memory_usage();
    const size_t database_memory_usage = index.get_database_memory_usage();
    const size_t ITLB = index.get_ITLB();
    const size_t index_ITLB = index.get_index_ITLB();
    const size_t num_nodes = index.get_num_nodes();
    const size_t num_inner_nodes = index.get_num_inner_nodes();
    const double num_values_per_leaf = index.get_num_values_per_leaf();

    reportfln("constr time in seconds: %d", constr_sec);
    reportfln("memory usage in MiB: %d", memory_usage / (1024.0 * 1024.0));
    reportfln("- index: %d", index_memory_usage / (1024.0 * 1024.0));
    reportfln("- database: %d", database_memory_usage / (1024.0 * 1024.0));
    reportfln("ITLB in MiB: %d", ITLB / (1024.0 * 1024.0));
    reportfln("- index: %d", index_ITLB / (1024.0 * 1024.0));
    reportfln("number of nodes: %d", num_nodes);
    reportfln("number of inner nodes: %d", num_inner_nodes);
    reportfln("number of values per leaf: %f", num_values_per_leaf);

    EXPERIMENT_APPEND("preprocess", {{"constr_sec", constr_sec},  //
                                     {"memory_usage", memory_usage},
                                     {"index_memory_usage", index_memory_usage},
                                     {"database_memory_usage", database_memory_usage},
                                     {"ITLB", ITLB},
                                     {"index_ITLB", index_ITLB},
                                     {"num_nodes", num_nodes},
                                     {"num_inner_nodes", num_inner_nodes},
                                     {"num_values_per_leaf", num_values_per_leaf}});

    /**
     *  Search queries
     */
    std::vector<ScoreVec> results(query_size);
    for (size_t i = 0; i < query_curves.size(); ++i) {
        results[i].reserve(1U << 10);
    }

    auto [min_range, max_range, step] = parse_range(hamming_ranges);

    for (uint32_t hamming_range = min_range; hamming_range <= max_range; hamming_range += step) {
        tfm::printfln("\nhamming range: %d", hamming_range);

        size_t num_candidates = 0;

        START_TIMER(searching);
        for (size_t i = 0; i < query_size; ++i) {
            results[i].clear();
            const Curve<Point>& query_curve = query_curves[i];
            num_candidates += index.range_search(query_curve, hamming_range, [&](size_t collision_pos) {
                results[i].push_back(std::make_pair(collision_pos, 1.0));
            });
        }
        STOP_TIMER_V(searching);

        const double ave_search_time_in_ms = GET_TIMER_MS(searching) / query_size;
        reportfln("search time per query in ms: %g", ave_search_time_in_ms);

#ifdef ENABLE_VERIFICATION
        size_t num_true_results = 0;

        START_TIMER(verifying);
        for (size_t i = 0; i < query_size; ++i) {
            const Curve<Point>& query_curve = query_curves[i];
            for (const Score& score : results[i]) {
                double frechet_dist = compute_frechet_discrete(query_curve, base_curves[score.first], frechet_range);
                if (0.0 <= frechet_dist && frechet_dist <= frechet_range) {
                    ++num_true_results;
                }
            }
        }
        STOP_TIMER_V(verifying);

        double ave_verify_time_in_ms = GET_TIMER_MS(verifying) / query_size;
        reportfln("verify time per query in ms: %g", ave_verify_time_in_ms);
#endif

        size_t num_results = 0;
        for (size_t i = 0; i < query_size; ++i) {
            num_results += results[i].size();
        }

        const double ave_results = double(num_results) / query_size;
        const double ave_candidates = double(num_candidates) / query_size;
        reportfln("number of results per query: %g (%g%%)", ave_results, ave_results / base_size * 100);
        reportfln("number of candidates (by multi-index) per query: %g (%g%%)",  //
                  ave_candidates, ave_candidates / base_size * 100);

#ifdef ENABLE_VERIFICATION
        const double ave_true_results = double(num_true_results) / query_size;
        reportfln("number of true results per query: %g (%g%%)", ave_true_results, ave_true_results / base_size * 100);

        EXPERIMENT_APPEND("search", {{"hamming_range", size_t(hamming_range)},  //
                                     {"ave_search_time_in_ms", ave_search_time_in_ms},
                                     {"ave_verify_time_in_ms", ave_verify_time_in_ms},
                                     {"ave_results", ave_results},
                                     {"ave_true_results", ave_true_results},
                                     {"ave_candidates", ave_candidates}});
#else
        EXPERIMENT_APPEND("search", {{"hamming_range", size_t(hamming_range)},  //
                                     {"ave_search_time_in_ms", ave_search_time_in_ms},
                                     {"ave_results", ave_results},
                                     {"ave_candidates", ave_candidates}});
#endif

        if (ABORT_THRESHOLD_MS < ave_search_time_in_ms) {
            warnfln("Stop search since the elapsed time exceeds ABORT_THRESHOLD_MS (%g)", ABORT_THRESHOLD_MS);
            break;
        }
    }

    tfm::printfln("saving %s", result_path);
    make_directory(result_dir);
    EXPERIMENT_SAVE_WITH_PATH(result_path);

    return 0;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    warnfln("The code is running in debug mode (or at least with assertions enabled)");
#endif

    cmdline::parser p;
    p.add<std::string>("base_path", 'b', "input file path of database", true);
    p.add<std::string>("query_path", 'q', "input file path of queries", true);
    p.add<std::string>("result_dir", 'o', "output directory path of results", false, "results");
    p.add<double>("frechet_range", 'r', "frechet range", true);
    p.add<std::string>("hamming_ranges", 'h', "hamming ranges (min:max:step)", false, "0:16:2");
    p.add<uint32_t>("repeatations", 'L', "repeatations of functions in LSH (i.e. code length)", false, 64);
    p.add<uint32_t>("buckets", 'B', "number of buckets", false, 8);
    p.add<uint32_t>("reduction_threshold", 'T', "threshold of node reduction", false, 1);
    p.add<double>("resolution_factor", 'F', "resolution factor in LSH", false, 4.0);
    p.add<uint64_t>("seed", 'S', "random seed in LSH", false, 114514);
    p.parse_check(argc, argv);

    auto base_path = p.get<std::string>("base_path");
    auto query_path = p.get<std::string>("query_path");
    size_t base_dimensions = read_dataset_dimensions(base_path);
    size_t query_dimensions = read_dataset_dimensions(query_path);

    if (base_dimensions != query_dimensions) {
        errorfln("base_dimensions != query_dimensions");
        return 1;
    }
    tfm::printfln("dimensionality %d", base_dimensions);

    const auto repeatations = p.get<uint32_t>("repeatations");

    if (base_dimensions == 1) {
        using point_type = Point1D;
        switch (repeatations) {
            case 32:
                return main_template<point_type, 32>(p);
            case 64:
                return main_template<point_type, 64>(p);
            default:
                errorfln("unsupported repeatations %d", repeatations);
        };
    } else if (base_dimensions == 2) {
        using point_type = Point2D;
        switch (repeatations) {
            case 32:
                return main_template<point_type, 32>(p);
            case 64:
                return main_template<point_type, 64>(p);
            default:
                errorfln("unsupported repeatations %d", repeatations);
        };
    } else {
        errorfln("unsupported dimensions %d", base_dimensions);
    }

    return 1;
}