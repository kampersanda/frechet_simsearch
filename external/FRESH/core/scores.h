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
#pragma once

#ifndef CURVEDIST_RECALL_H
#define CURVEDIST_RECALL_H

#include <boost/filesystem.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "io.h"
#include "msgpack.hpp"
#include "prelude.h"
#include "types.h"

typedef std::pair<size_t, size_t> id_pair_t;

struct Baseline {
    Baseline() {}

    std::string queries_path;
    std::string data_path;
    double range;
    std::set<std::pair<size_t, size_t>> pairs;

    MSGPACK_DEFINE_MAP(queries_path, data_path, range, pairs);
};

std::string baseline_path(const std::string basedir, std::string queries_path, std::string data_path, double range) {
    std::replace(queries_path.begin(), queries_path.end(), '/', '_');
    std::replace(queries_path.begin(), queries_path.end(), '.', '_');
    std::replace(data_path.begin(), data_path.end(), '/', '_');
    std::replace(data_path.begin(), data_path.end(), '.', '_');
    std::stringstream strstr;
    strstr << basedir << "/"
           << "baseline"
           << "-" << queries_path << "-" << data_path << "-" << range << ".mpk";
    return strstr.str();
}

void dump_baseline(const std::set<std::pair<size_t, size_t>>& pairs, const std::string basedir,
                   const std::string queries_path, const std::string data_path, double range) {
    Baseline b;
    b.pairs = pairs;
    b.queries_path = queries_path;
    b.data_path = data_path;
    b.range = range;

    std::string path = baseline_path(basedir, queries_path, data_path, range);
    std::cout << "Output path " << path << std::endl;
    if (boost::filesystem::exists(path)) {
        std::cout << "Baseline already exists, skip dumping" << std::endl;
        return;
    }
    boost::filesystem::path base_path(basedir);
    if (!boost::filesystem::exists(base_path)) {
        if (!boost::filesystem::create_directory(base_path)) {
            throw std::logic_error("Unable to create output directory");
        }
    }

    std::ofstream out(path);
    msgpack::pack(out, b);
}

Baseline load_baseline(const std::string basedir, const std::string queries_path, const std::string data_path,
                       double range) {
    std::string path = baseline_path(basedir, queries_path, data_path, range);
    std::cout << "loading baseline from " << path << std::endl;

    std::ifstream input_stream(path);
    if (!input_stream.good()) {
        throw std::logic_error("Could not open file");
    }

    // read the file in memory
    std::stringstream strstr;
    strstr << input_stream.rdbuf();

    std::string buffer = strstr.str();
    msgpack::object_handle oh = msgpack::unpack(buffer.data(), buffer.size());
    msgpack::object obj = oh.get();
    Baseline b;
    obj.convert(b);
    return b;
}

double recall(const std::set<std::pair<size_t, size_t>>& predictions,
              const std::set<std::pair<size_t, size_t>>& baseline) {
    size_t true_positives = 0;
    for (const auto& positive : baseline) {
        if (predictions.count(positive) > 0) {
            true_positives++;
        }
    }
    printf("Counted true positives %lu\n", true_positives);
    return true_positives / ((double)baseline.size());
}

double precision(const std::set<std::pair<size_t, size_t>>& predictions,
                 const std::set<std::pair<size_t, size_t>>& baseline) {
    size_t true_positives = 0;
    for (const auto& positive : predictions) {
        if (baseline.count(positive) > 0) {
            true_positives++;
        }
    }
    printf("Predicted positives %lu\n", predictions.size());
    return true_positives / ((double)predictions.size());
}

double recall(const std::vector<std::pair<double, std::pair<size_t, size_t>>>& predictions,
              const std::set<std::pair<size_t, size_t>>& baseline) {
    size_t true_positives = 0;
    for (const auto& prediction : predictions) {
        if (prediction.first > 0 && baseline.count(prediction.second) > 0) {
            true_positives++;
        }
    }
    printf("Counted true positives %lu\n", true_positives);
    return true_positives / ((double)baseline.size());
}

double precision(const std::vector<std::pair<double, std::pair<size_t, size_t>>>& predictions,
                 const std::set<std::pair<size_t, size_t>>& baseline) {
    size_t true_positives = 0;
    size_t predicted_positives = 0;
    for (const auto& prediction : predictions) {
        if (prediction.first > 0) {
            predicted_positives++;
            if (baseline.count(prediction.second) > 0) {
                true_positives++;
            }
        }
    }
    printf("Predicted positives %lu (of which true: %lu)\n", predicted_positives, true_positives);
    return true_positives / ((double)predicted_positives);
}

template <typename Point>
void print_false_positives(const std::set<std::pair<size_t, size_t>>& predictions,
                           const std::set<std::pair<size_t, size_t>>& baseline,
                           const std::vector<Curve<Point>>& queries, const std::vector<Curve<Point>>& dataset) {
    std::vector<double> tmp1, tmp2;

    std::vector<std::pair<double, std::pair<size_t, size_t>>> false_positives;
    std::cout << "== false positives ==" << std::endl;
    for (const auto& positive : predictions) {
        if (baseline.count(positive) == 0) {
            auto query = std::find_if(queries.begin(), queries.end(),
                                      [&](const Curve<Point>& c) { return c.id == positive.first; });
            auto curve = std::find_if(dataset.begin(), dataset.end(),
                                      [&](const Curve<Point>& c) { return c.id == positive.second; });
            double d = discrete_frechet_distance(query->points, curve->points, tmp1, tmp2);
            false_positives.emplace_back(d, std::make_pair(query->id, curve->id));
        }
    }
    std::sort(false_positives.begin(), false_positives.end());
    for (auto& fp : false_positives) {
        std::cout << fp.first << " (" << fp.second.first << ", " << fp.second.second << ")" << std::endl;
    }
    std::cout << "====" << std::endl;
}

template <typename Point>
std::set<std::pair<size_t, size_t>> result_as_pairs(const std::vector<std::vector<size_t>>& results,
                                                    const std::vector<Curve<Point>>& queries,
                                                    const std::vector<Curve<Point>>& dataset) {
    if (results.size() != queries.size()) {
        throw std::logic_error("Result of the wrong size");
    }
    std::set<std::pair<size_t, size_t>> pairs;
    for (size_t query_idx = 0; query_idx < results.size(); query_idx++) {
        size_t query_curve_id = queries.at(query_idx).id;
        for (const size_t data_idx : results.at(query_idx)) {
            size_t data_curve_id = dataset.at(data_idx).id;
            pairs.insert(std::make_pair(query_curve_id, data_curve_id));
        }
    }

    return pairs;
};

typedef std::pair<double, double> double_pair_t;

std::vector<std::pair<double, double>> precision_recall_curve(
    const std::vector<std::pair<double, id_pair_t>>& scored_results, const std::set<id_pair_t>& baseline,
    const size_t total_pairs) {
    // flags for the pairs
    std::vector<bool> flags(scored_results.size());

    // find the distinct scores
    std::set<double> scores_set;
    for (size_t i = 0; i < scored_results.size(); i++) {
        auto elem = scored_results[i];
        scores_set.insert(elem.first);
        if (baseline.count(elem.second) > 0) {
            flags[i] = true;
        } else {
            flags[i] = false;
        }
    }

    // build the curve
    std::vector<std::pair<double, double>> pr_curve;
    pr_curve.reserve(scores_set.size() + 1);
    pr_curve.emplace_back(0, 1);
    pr_curve.emplace_back(1, 0);
    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    printf("total pairs %zu\n", total_pairs);
    // printf("total pairs %d\n", total_pairs);
    const size_t num_positives = baseline.size();
    const size_t num_negatives = total_pairs - num_positives;
    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    printf("positives %zu negative %zu (solution pairs %zu)\n", num_positives, num_negatives, scored_results.size());
    // printf("positives %d negative %d (solution pairs %d)\n", num_positives, num_negatives, scored_results.size());
    for (double thresh : scores_set) {
        size_t tp = 0;
        size_t fp = 0;
        for (size_t i = 0; i < scored_results.size(); i++) {
            auto elem = scored_results[i];
            if (elem.first >= thresh) {
                if (flags[i]) {
                    tp++;
                } else {
                    fp++;
                }
            }
        }
        double recall = tp / ((double)num_positives);
        double precision = tp / ((double)tp + fp);
        pr_curve.emplace_back(recall, precision);
    }
    std::sort(pr_curve.begin(), pr_curve.end(), [&](double_pair_t a, double_pair_t b) {
        if (a.first < b.first) {
            return true;
        } else if (a.first == b.first) {
            return a.second > b.second;  // We want them decreasing in the second element
        }
        return false;
    });

    return pr_curve;
}

std::vector<std::pair<double, double>> roc(const std::vector<std::pair<double, id_pair_t>>& scored_results,
                                           const std::set<id_pair_t>& baseline, const size_t total_pairs) {
    // flags for the pairs
    std::vector<bool> flags(scored_results.size());

    // find the distinct scores
    std::set<double> scores_set;
    /* std::vector<double> scores_list; */
    /* scores_list.reserve(scored_results.size()); */
    for (size_t i = 0; i < scored_results.size(); i++) {
        auto elem = scored_results[i];
        scores_set.insert(elem.first);
        /* if (elem.first > 0) { */
        /*   scores_list.push_back(elem.first); */
        /* } */
        if (baseline.count(elem.second) > 0) {
            flags[i] = true;
        } else {
            flags[i] = false;
        }
    }
    /* std::sort(scores_list.begin(), scores_list.end()); */
    /* EXPERIMENT_APPEND("scores_stats", { */
    /*     {"min", scores_list[0]}, */
    /*     {"max", scores_list[scores_list.size() -1]}, */
    /*     {"median", scores_list[scores_list.size()/2]}, */
    /*     {"first quartile", scores_list[scores_list.size()/4]}, */
    /*     {"third quartile", scores_list[3*scores_list.size()/4]} */
    /*     }); */

    // build the ROC
    std::vector<std::pair<double, double>> roc_curve;
    roc_curve.reserve(scores_set.size() + 1);
    roc_curve.emplace_back(0, 0);
    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    printf("total pairs %zu\n", total_pairs);
    // printf("total pairs %d\n", total_pairs);
    const size_t num_positives = baseline.size();
    const size_t num_negatives = total_pairs - num_positives;
    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    printf("positives %zu negative %zu (solution pairs %zu)\n", num_positives, num_negatives, scored_results.size());
    // printf("positives %d negative %d (solution pairs %d)\n", num_positives, num_negatives, scored_results.size());
    for (double thresh : scores_set) {
        size_t tp = 0;
        size_t fp = 0;
        for (size_t i = 0; i < scored_results.size(); i++) {
            auto elem = scored_results[i];
            if (elem.first >= thresh) {
                if (flags[i]) {
                    tp++;
                } else {
                    fp++;
                }
            }
        }
        double tpr = tp / ((double)num_positives);
        double fpr = fp / ((double)num_negatives);
        /* printf("fp=%d (%f)   tp=%d (%f)\n", fp, tp, fpr, tpr); */
        roc_curve.emplace_back(fpr, tpr);
    }
    std::sort(roc_curve.begin(), roc_curve.end());

    return roc_curve;
}

void report_scores(const std::vector<std::pair<double, id_pair_t>>& scored_results,
                   const std::set<id_pair_t>& baseline) {
    size_t true_positives_cnt = 0, false_positives_cnt = 0;
    std::map<double, size_t> true_positives, false_positives;
    for (auto prediction : scored_results) {
        if (prediction.first > 0) {
            if (baseline.count(prediction.second) > 0) {
                true_positives_cnt++;
                if (true_positives.count(prediction.first) > 0) {
                    true_positives[prediction.first] += 1;
                } else {
                    true_positives[prediction.first] = 1;
                }
            } else {
                false_positives_cnt++;
                if (false_positives.count(prediction.first) > 0) {
                    false_positives[prediction.first] += 1;
                } else {
                    false_positives[prediction.first] = 1;
                }
            }
        }
    }
    /* printf("True positives scores\n"); */
    for (auto e : true_positives) {
        double frac = e.second / ((double)true_positives_cnt);
        /* std::cout << e.first << " " << e.second << " (" << frac << ")" << std::endl; */
        EXPERIMENT_APPEND("true_positives_scores", {{"score", e.first}, {"count", e.second}, {"fraction", frac}});
    }
    /* printf("False positives scores\n"); */
    for (auto e : false_positives) {
        double frac = e.second / ((double)false_positives_cnt);
        /* std::cout << e.first << " " << e.second << " (" << frac << ")" << std::endl; */
        EXPERIMENT_APPEND("false_positives_scores", {{"score", e.first}, {"count", e.second}, {"fraction", frac}});
    }
}

double auc(const std::vector<std::pair<double, double>>& roc_curve) {
    double area = 0.0;
    for (size_t i = 0; i < roc_curve.size() - 1; i++) {
        auto a = roc_curve[i];
        auto b = roc_curve[i + 1];
        double h = b.first - a.first;
        double bB = b.second + a.second;
        area += h * bB / 2.0;
    }
    return area;
}

double auc(const double fpr, const double tpr) {
    std::vector<std::pair<double, double>> roc_curve;
    roc_curve.emplace_back(0, 0);
    roc_curve.emplace_back(fpr, tpr);
    roc_curve.emplace_back(1, 1);
    return auc(roc_curve);
}

void dump_result(const std::vector<std::pair<double, std::pair<size_t, size_t>>>& result,
                 const std::set<std::pair<size_t, size_t>>& baseline, const std::string& output_path) {
    std::vector<std::pair<size_t, size_t>> false_negatives;
    std::vector<std::pair<size_t, size_t>> true_positives;
    std::vector<std::pair<size_t, size_t>> false_positives;
    std::set<std::pair<size_t, size_t>> prediction_set;
    size_t false_negatives_cnt = 0;
    for (auto& pred : result) {
        if (pred.first > 0) {
            prediction_set.insert(pred.second);
            if (baseline.count(pred.second) == 0) {
                false_positives.push_back(pred.second);
            }
        }
    }

    for (auto& positive : baseline) {
        if (prediction_set.count(positive) == 0) {
            false_negatives.push_back(positive);
            false_negatives_cnt++;
        } else {
            true_positives.push_back(positive);
        }
    }
    std::cout << "False negatives " << false_negatives_cnt << std::endl;

    std::ofstream out_stream(output_path);
    msgpack::packer<std::ofstream> pk(out_stream);
    pk.pack_map(3);
    pk.pack(std::string("false_negatives"));
    pk.pack(false_negatives);
    pk.pack(std::string("true_positives"));
    pk.pack(true_positives);
    pk.pack(std::string("false_positives"));
    pk.pack(false_positives);
}

#endif  // CURVEDIST_RECALL_H
