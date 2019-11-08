#include <cmath>
#include <ctime>
#include <random>

#include <FRESH/core/experiment_reporter.h>
#include <FRESH/core/scores.h>
#include <FRESH/core/timer.h>

#include <cmdline.h>

#include <frechet_tools.hpp>
#include <fresh_index.hpp>
#include <groundtruth.hpp>

using namespace frechet_simsearch;

std::string make_result_path(std::string result_dir, std::string base_path, std::string query_path,
                             double frechet_range, std::string hamming_similarities, Configure cfg) {
    auto result_path = tfm::format("%s/fresh_score-%s-%s-%g-%s-%s.mpk",  //
                                   result_dir, normalize_filepath(base_path), normalize_filepath(query_path),  //
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
    const auto groundtruth_dir = p.get<std::string>("groundtruth_dir");
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
     *  Load groundtruth
     */
    const auto groundtruth_path = get_range_groundtruth_path(groundtruth_dir, query_path, base_path, frechet_range);
    if (!boost::filesystem::exists(groundtruth_path)) {
        errorfln("groundtruth data does not exist: %s", groundtruth_path);
        return 1;
    }
    const RangeGroundTruth gt = load_range_groundtruth(groundtruth_dir, query_path, base_path, frechet_range);

    /**
     *  Build Index
     */
    START_TIMER(preprocessing);
    FreshIndex<Point> index(base_curves, cfg);
    STOP_TIMER_V(preprocessing);

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

        size_t num_candidates = 0;

        START_TIMER(searching);
        for (size_t i = 0; i < query_size; ++i) {
            results[i].clear();

            const Curve<Point>& query_curve = query_curves[i];
            index.count_collisions(query_curve, counters);

            for (size_t collision_pos = 0; collision_pos < base_size; ++collision_pos) {
                if (counters[collision_pos] > similarity) {
                    const Curve<Point>& base_curve = base_curves[collision_pos];
                    const double frechet_dist = compute_frechet_discrete(query_curve, base_curve, frechet_range);
                    if (0.0 <= frechet_dist && frechet_dist <= frechet_range) {
                        results[i].push_back(std::make_pair(collision_pos, frechet_dist));
                    }
                    ++num_candidates;
                }
            }
        }
        STOP_TIMER_V(searching);

        size_t num_results = 0;
        for (size_t i = 0; i < query_size; ++i) {
            num_results += results[i].size();
        }

        const double ave_results = double(num_results) / query_size;
        const double ave_candidates = double(num_candidates) / query_size;

        tfm::printfln("number of results per query: %g (%g%%)", ave_results, ave_results / base_size * 100);
        tfm::printfln("number of candidates per query: %g (%g%%)", ave_candidates, ave_candidates / base_size * 100);
        tfm::printfln("candidates / results = %g", ave_candidates / ave_results);

        auto score_map = make_score_map(results, query_curves, base_curves);
        const double recall_score = compute_recall(score_map, gt.score_map);
        const double precision_score = num_results / double(num_candidates);
        const double fscore = 2 * precision_score * recall_score / (precision_score + recall_score);

        reportfln("recall score: %g", recall_score);  // false negative ratio
        reportfln("precision score: %g", precision_score);
        reportfln("F-score: %g", fscore);

        EXPERIMENT_APPEND("score", {{"hamming_similarity", size_t(similarity)},  //
                                    {"ave_results", ave_results},  //
                                    {"ave_candidates", ave_candidates},  //
                                    {"recall_score", recall_score},  //
                                    {"precision_score", precision_score},
                                    {"fscore", fscore}});
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
    p.add<std::string>("groundtruth_dir", 'g', "input directory path of true solutions", false, "groundtruth");
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