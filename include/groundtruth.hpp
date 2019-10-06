#pragma once

#include <FRESH/core/scores.h>

#include "basic.hpp"

using Score = std::pair<size_t, double>;  // base_id, score
using ScoreVec = std::vector<Score>;
using ScoreMap = std::map<size_t, ScoreVec>;  // query_id -> (base_id, score)s

struct RangeGroundTruth {
    RangeGroundTruth() = default;

    std::string query_path;
    std::string base_path;
    double range;
    ScoreMap score_map;  // query_id -> (base_id, score)s

    MSGPACK_DEFINE_MAP(query_path, base_path, range, score_map);
};

struct TopkGroundTruth {
    TopkGroundTruth() = default;

    std::string query_path;
    std::string base_path;
    size_t topk;
    ScoreMap score_map;  // query_id -> (base_id, score)s

    MSGPACK_DEFINE_MAP(query_path, base_path, topk, score_map);
};

inline std::string get_range_groundtruth_path(std::string groundtruth_dir,  //
                                              std::string query_path, std::string base_path, double range) {
    query_path = normalize_filepath(query_path);
    base_path = normalize_filepath(base_path);
    auto range_str = normalize_filepath(tfm::format("%g", range));
    return tfm::format("%s/range_groundtruth-%s-%s-%s.mpk", groundtruth_dir, query_path, base_path, range_str);
}

inline std::string get_topk_groundtruth_path(std::string groundtruth_dir,  //
                                             std::string query_path, std::string base_path, size_t topk) {
    query_path = normalize_filepath(query_path);
    base_path = normalize_filepath(base_path);
    return tfm::format("%s/topk_groundtruth-%s-%s-%d.mpk", groundtruth_dir, query_path, base_path, topk);
}

inline void dump_range_groundtruth(const ScoreMap& score_map, std::string groundtruth_dir,  //
                                   std::string query_path, std::string base_path, double range) {
    RangeGroundTruth gt;
    gt.query_path = query_path;
    gt.base_path = base_path;
    gt.range = range;
    gt.score_map = score_map;

    auto groundtruth_path = get_range_groundtruth_path(groundtruth_dir, query_path, base_path, range);
    tfm::printfln("output groundtruth %s", groundtruth_path);

    if (boost::filesystem::exists(groundtruth_path)) {
        warnfln("groundtruth data already exists, so skip dumping");
        return;
    }
    make_directory(groundtruth_dir);

    std::ofstream out(groundtruth_path);
    msgpack::pack(out, gt);
}

inline void dump_topk_groundtruth(const ScoreMap& score_map, std::string groundtruth_dir,  //
                                  std::string query_path, std::string base_path, size_t topk) {
    TopkGroundTruth gt;
    gt.query_path = query_path;
    gt.base_path = base_path;
    gt.topk = topk;
    gt.score_map = score_map;

    auto groundtruth_path = get_topk_groundtruth_path(groundtruth_dir, query_path, base_path, topk);
    tfm::printfln("output groundtruth %s", groundtruth_path);

    if (boost::filesystem::exists(groundtruth_path)) {
        warnfln("groundtruth data already exists, so skip dumping");
        return;
    }
    make_directory(groundtruth_dir);

    std::ofstream out(groundtruth_path);
    msgpack::pack(out, gt);
}

inline RangeGroundTruth  //
load_range_groundtruth(std::string groundtruth_dir, std::string query_path, std::string base_path, double range) {
    auto groundtruth_path = get_range_groundtruth_path(groundtruth_dir, query_path, base_path, range);
    tfm::printfln("load groundtruth %s", groundtruth_path);

    std::ifstream input_stream(groundtruth_path);
    if (!input_stream.good()) {
        errorfln("could not open %s", groundtruth_path);
        exit(1);
    }

    // read the file in memory
    std::stringstream strstr;
    strstr << input_stream.rdbuf();

    std::string buffer = strstr.str();
    msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());
    msgpack::object obj = oh.get();

    RangeGroundTruth gt;
    obj.convert(gt);
    return gt;
}

inline TopkGroundTruth  //
load_topk_groundtruth(std::string groundtruth_dir, std::string query_path, std::string base_path, size_t topk) {
    auto groundtruth_path = get_topk_groundtruth_path(groundtruth_dir, query_path, base_path, topk);
    tfm::printfln("load groundtruth %s", groundtruth_path);

    std::ifstream input_stream(groundtruth_path);
    if (!input_stream.good()) {
        errorfln("could not open %s", groundtruth_path);
        exit(1);
    }

    // read the file in memory
    std::stringstream strstr;
    strstr << input_stream.rdbuf();

    std::string buffer = strstr.str();
    msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());
    msgpack::object obj = oh.get();

    TopkGroundTruth gt;
    obj.convert(gt);
    return gt;
}

template <class Point>
inline ScoreMap make_score_map(const std::vector<ScoreVec>& results,  //
                               const std::vector<Curve<Point>>& query_curves,
                               const std::vector<Curve<Point>>& base_curves) {
    if (results.size() != query_curves.size()) {
        errorfln("ScoreVec of the wrong size");
        exit(1);
    }

    ScoreVec scores;
    ScoreMap score_map;

    for (size_t query_pos = 0; query_pos < query_curves.size(); ++query_pos) {
        size_t query_id = query_curves[query_pos].id;

        scores.clear();
        for (const Score& score : results[query_pos]) {
            size_t base_id = base_curves[score.first].id;
            scores.push_back(std::make_pair(base_id, score.second));
        }
        score_map.insert(std::make_pair(query_id, scores));
    }

    return score_map;
}

inline double compute_recall(const ScoreMap& predictions, const ScoreMap& groundtruth) {
    size_t true_positives = 0;
    size_t true_solutions = 0;
    for (const auto& positive : groundtruth) {
        const auto& pred_itr = predictions.find(positive.first);
        if (pred_itr == predictions.end()) {
            errorfln("query id in groundtruth is not contained in predictions");
            exit(1);
        }
        const ScoreVec& pred_scores = pred_itr->second;
        for (std::pair<size_t, double> score : positive.second) {
            size_t base_id = score.first;
            auto it = std::find_if(pred_scores.begin(), pred_scores.end(),
                                   [base_id](const std::pair<size_t, double>& x) { return x.first == base_id; });
            if (it != pred_scores.end()) {
                true_positives++;
            }
        }
        true_solutions += positive.second.size();
    }
    return true_positives / double(true_solutions);
}

inline double compute_topk_recall(const ScoreMap& predictions, const ScoreMap& groundtruth, size_t topk) {
    size_t true_positives = 0;
    size_t true_solutions = 0;
    for (const auto& positive : groundtruth) {
        const auto& pred_itr = predictions.find(positive.first);
        if (pred_itr == predictions.end()) {
            errorfln("query id in groundtruth is not contained in predictions");
            exit(1);
        }

        const ScoreVec& pred_scores = pred_itr->second;

        if (pred_scores.size() < topk) {
            warnfln("number of predictions expected is %d but actually %d", topk, pred_scores.size());
        }

        const auto pred_scores_end = pred_scores.begin() + std::min(topk, pred_scores.size());

        for (std::pair<size_t, double> score : positive.second) {
            size_t base_id = score.first;
            auto it = std::find_if(pred_scores.begin(), pred_scores_end,
                                   [base_id](const std::pair<size_t, double>& x) { return x.first == base_id; });
            if (it != pred_scores_end) {
                true_positives++;
            }
        }
        true_solutions += positive.second.size();
    }
    return true_positives / double(true_solutions);
}

inline std::vector<std::pair<size_t, size_t>>  //
get_true_positives(const ScoreMap& predictions, const ScoreMap& groundtruth) {
    std::vector<std::pair<size_t, size_t>> true_positives;
    for (const auto& positive : groundtruth) {
        const auto& pred_itr = predictions.find(positive.first);
        if (pred_itr == predictions.end()) {
            errorfln("query id in groundtruth is not contained in predictions");
            exit(1);
        }
        const ScoreVec& pred_scores = pred_itr->second;
        for (std::pair<size_t, double> score : positive.second) {
            size_t base_id = score.first;
            auto it = std::find_if(pred_scores.begin(), pred_scores.end(),
                                   [base_id](const std::pair<size_t, double>& x) { return x.first == base_id; });
            if (it != pred_scores.end()) {
                true_positives.push_back(std::make_pair(base_id, it->first));
            }
        }
    }
    return true_positives;
}