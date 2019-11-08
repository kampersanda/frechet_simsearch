#include <cmath>
#include <ctime>
#include <random>

#include <FRESH/core/experiment_reporter.h>
#include <FRESH/core/scores.h>
#include <FRESH/core/timer.h>

#include <cmdline.h>

#include <frechet_hasher.hpp>
#include <groundtruth.hpp>

using namespace frechet_simsearch;

extern std::vector<ScoreVec> g_results;
#pragma omp threadprivate(g_results)
std::vector<ScoreVec> g_results;

template <class Point>
int main_template(const cmdline::parser& p) {
    tfm::printfln("num threads: %d", omp_get_max_threads());

    const auto base_path = p.get<std::string>("base_path");
    const auto query_path = p.get<std::string>("query_path");
    const auto groundtruth_dir = p.get<std::string>("groundtruth_dir");
    const auto frechet_range = p.get<double>("frechet_range");

    const std::vector<Curve<Point>> base_curves = load_dataset_msgpack<Point>(base_path);
    const std::vector<Curve<Point>> query_curves = load_dataset_msgpack<Point>(query_path);
    tfm::printfln("num base curves : %d", base_curves.size());
    tfm::printfln("num query curves: %d", query_curves.size());

#pragma omp parallel
    { g_results.resize(query_curves.size()); }

    size_t progress = 0;

    START_TIMER(searching);
#pragma omp parallel for
    for (size_t query_pos = 0; query_pos < query_curves.size(); ++query_pos) {
        const Curve<Point>& query_curve = query_curves[query_pos];
        for (size_t base_pos = 0; base_pos < base_curves.size(); ++base_pos) {
            const Curve<Point>& base_curve = base_curves[base_pos];
            double distance = compute_frechet_discrete(query_curve, base_curve, frechet_range);
            if (0.0 <= distance && distance <= frechet_range) {
                g_results[query_pos].push_back(std::make_pair(base_pos, distance));
            }
        }
#pragma omp critical
        {
            ++progress;
            if (progress % 100 == 0) {
                tfm::printfln("%d/%d", progress, query_curves.size());
            }
        }
    }
    STOP_TIMER_V(searching);

    const double search_ms_per_query = GET_TIMER_MS(searching) / query_curves.size();
    reportfln("search time in milliseconds per query: %g", search_ms_per_query);

    std::vector<ScoreVec> results(query_curves.size());
#pragma omp parallel for
    for (size_t i = 0; i < query_curves.size(); ++i) {
        std::copy(g_results[i].begin(), g_results[i].end(), std::back_inserter(results[i]));
    }

    size_t num_results = 0;
    for (const ScoreVec& res : results) {
        num_results += res.size();
    }
    const double ave_results = double(num_results) / results.size();
    reportfln("number of results per query: %g (%g%%)", ave_results, ave_results / base_curves.size() * 100);

    std::string groundtruth_path = get_range_groundtruth_path(groundtruth_dir, query_path, base_path, frechet_range);

    if (!boost::filesystem::exists(groundtruth_path)) {  // need to generate?
        tfm::printfln("Generating true solutions...");
        auto score_map = make_score_map(results, query_curves, base_curves);
        dump_range_groundtruth(score_map, groundtruth_dir, query_path, base_path, frechet_range);
    } else {
        warnfln("groundtruth data already exists, so skip dumping");
    }

    return 0;
}

int main(int argc, char** argv) {
#ifndef NDEBUG
    warnfln("The code is running in debug mode (or at least with assertions enabled)");
#endif

    cmdline::parser p;
    p.add<std::string>("base_path", 'b', "input file path of database", true);
    p.add<std::string>("query_path", 'q', "input file path of queries", true);
    p.add<std::string>("groundtruth_dir", 'g', "output directory path of true solutions", false, "groundtruth");
    p.add<double>("frechet_range", 'r', "frechet range", true);
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