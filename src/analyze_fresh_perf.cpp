#include <cmath>
#include <ctime>
#include <random>

#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include <FRESH/core/experiment_reporter.h>
#include <FRESH/core/scores.h>
#include <FRESH/core/timer.h>

#include <cmdline.h>

#include <frechet_tools.hpp>
#include <fresh_index.hpp>
#include <groundtruth.hpp>

#define ENABLE_VERIFICATION

// From Cedar (http://www.tkl.iis.u-tokyo.ac.jp/~ynaga/cedar/)
inline size_t get_process_size() {
#ifdef __APPLE__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&t_info), &t_info_count);
    return t_info.resident_size;
#else
    FILE* fp = std::fopen("/proc/self/statm", "r");
    size_t dummy(0), vm(0);
    std::fscanf(fp, "%ld %ld ", &dummy, &vm);  // get resident (see procfs)
    std::fclose(fp);
    return vm * ::getpagesize();
#endif
}

std::string make_result_path(std::string result_dir, std::string base_path, std::string query_path,
                             double frechet_range, std::string hamming_similarities, Configure cfg) {
    auto result_path = tfm::format("%s/fresh_perf-%s-%s-%g-%s-%s.mpk",  //
                                   result_dir, normalize_filepath(base_path), normalize_filepath(query_path),
                                   frechet_range, normalize_filepath(hamming_similarities), cfg.to_str());
#ifdef EXPERIMENT_COMPRESS
    result_path += ".bz2";
#endif
    return result_path;
}

template <class Point>
int main_template(const cmdline::parser& p) {
    tfm::printfln("num threads: %d", omp_get_max_threads());

    const auto base_path = p.get<std::string>("base_path");
    const auto query_path = p.get<std::string>("query_path");
    const auto result_dir = p.get<std::string>("result_dir");
    const auto frechet_range = p.get<double>("frechet_range");
    const auto hamming_similarities = p.get<std::string>("hamming_similarities");
    const auto concatenations = p.get<uint32_t>("concatenations");
    const auto repeatations = p.get<uint32_t>("repeatations");
    const auto resolution_factor = p.get<double>("resolution_factor");
    const auto seed = p.get<uint64_t>("seed");

    const double resolution = resolution_factor * frechet_range * Point::dimensions();

    tfm::printfln("frechet_range: %g", frechet_range);
    tfm::printfln("resolution: %g (factor: %g)", resolution, resolution_factor);

    Configure cfg;
    {
        cfg.repeatations = repeatations;
        cfg.concatenations = concatenations;
        cfg.resolution = resolution;
        cfg.resolution_factor = resolution_factor;
        cfg.seed = seed;
    }
    cfg.show(std::cout);

    const std::string result_path = make_result_path(result_dir, base_path, query_path,  //
                                                     frechet_range, hamming_similarities, cfg);
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

    EXPERIMENT_TAG("algorithm", std::string("FRESH"));
    EXPERIMENT_TAG("base_path", base_path);
    EXPERIMENT_TAG("query_path", query_path);
    EXPERIMENT_TAG("base_size", base_size);
    EXPERIMENT_TAG("query_size", query_size);
    EXPERIMENT_TAG("frechet_range", frechet_range);
    EXPERIMENT_TAG("hamming_similarities", hamming_similarities);
    EXPERIMENT_TAG("concatenations", size_t(concatenations));
    EXPERIMENT_TAG("repeatations", size_t(repeatations));
    EXPERIMENT_TAG("resolution_factor", resolution_factor);
    EXPERIMENT_TAG("resolution", resolution);
    EXPERIMENT_TAG("seed", size_t(seed));
    EXPERIMENT_TAG("dimensions", Point::dimensions());

    /**
     *  Build index
     */
    START_TIMER(preprocessing);
    const size_t initial_process_size = get_process_size();
    FreshIndex<Point> index(base_curves, cfg);
    const size_t memory_usage = get_process_size() - initial_process_size;
    STOP_TIMER_V(preprocessing);

    const double constr_sec = GET_TIMER_MS(preprocessing) / 1000.0;

    reportfln("constr time in seconds: %d", constr_sec);
    reportfln("memory usage in MiB: %d", memory_usage / (1024.0 * 1024.0));

    EXPERIMENT_APPEND("preprocess", {{"constr_sec", constr_sec},  //
                                     {"memory_usage", memory_usage}});

    /**
     *  Search queries
     */
    std::vector<size_t> counters(base_size);
    std::vector<ScoreVec> results(query_size);
    for (size_t i = 0; i < query_size; ++i) {
        results[i].reserve(1U << 10);
    }

    auto [min_sim, max_sim, step] = parse_range(hamming_similarities);

    for (uint32_t similarity = min_sim; similarity <= max_sim; similarity += step) {
        tfm::printfln("\nhamming similarity: %d", similarity);

        START_TIMER(searching);
        for (size_t i = 0; i < query_size; ++i) {
            results[i].clear();

            const Curve<Point>& query_curve = query_curves[i];
            index.count_collisions(query_curve, counters);

            for (size_t collision_pos = 0; collision_pos < base_size; ++collision_pos) {
                if (counters[collision_pos] > similarity) {
                    results[i].push_back(std::make_pair(collision_pos, 1.0));
                }
            }
        }
        STOP_TIMER_V(searching);

        double ave_search_time_in_ms = GET_TIMER_MS(searching) / query_size;
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
        reportfln("number of results per query: %g (%g%%)", ave_results, ave_results / base_size * 100);

#ifdef ENABLE_VERIFICATION
        const double ave_true_results = double(num_true_results) / query_size;
        reportfln("number of true results per query: %g (%g%%)", ave_true_results, ave_true_results / base_size * 100);

        EXPERIMENT_APPEND("search", {{"hamming_similarity", size_t(similarity)},  //
                                     {"ave_search_time_in_ms", ave_search_time_in_ms},
                                     {"ave_verify_time_in_ms", ave_verify_time_in_ms},
                                     {"ave_results", ave_results},
                                     {"ave_true_results", ave_true_results}});
#else
        EXPERIMENT_APPEND("search", {{"hamming_similarity", size_t(similarity)},  //
                                     {"ave_search_time_in_ms", ave_search_time_in_ms},
                                     {"ave_results", ave_results}});
#endif
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
    p.add<std::string>("hamming_similarities", 'h', "hamming similarities (exclusive) (min:max:step)", false, "0:32:4");
    p.add<uint32_t>("concatenations", 'K', "concatenations of functions in LSH", false, 1);
    p.add<uint32_t>("repeatations", 'L', "repeatations of functions in LSH", false, 64);
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

    if (base_dimensions == 1) {
        return main_template<Point1D>(p);
    } else if (base_dimensions == 2) {
        return main_template<Point2D>(p);
    } else {
        errorfln("unsupported dimensions %d", base_dimensions);
    }

    return 1;
}