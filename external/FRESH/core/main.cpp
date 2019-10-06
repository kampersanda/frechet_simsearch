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
#include "frechet.h"
#include "hash.h"
#include "io.h"
#include "prelude.h"
#include "query.h"
#include "experiment_reporter.h"
#include "stats.h"
#include "timer.h"
#include "git_info.h"
#include "scores.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "rand.h"
#include <random>
#include <functional>

int parse_opts(int argc, char **argv,
               boost::program_options::variables_map &vm) {
  namespace po = boost::program_options;

  po::options_description desc("Experiment runner for range queries");
  desc.add_options()
          ("help", "print help message")
          ("algorithm", po::value<std::string>()->required(), "the algorithm to use for answering queries")
          ("dataset", po::value<std::string>(), "the dataset to be used")
          ("queries", po::value<std::string>(), "the queries to be answered")
          ("baseline", po::value<std::string>()->default_value("baselines"), "the queries to be answered")
          ("range", po::value<double>()->required(), "the range of the queries")
          ("seed", po::value<uint64_t>(), "the seed for the random generator")
          ("k", po::value<size_t>()->default_value(1))
          ("L", po::value<size_t>()->default_value(0))
          ("resolution-factor", po::value<double>())
          ("eval-dist", po::value<std::string>()->default_value("yes"), "evaluate the distance")
          ("dump-solution", po::value<std::string>(), "Dump the solution in msgpack format")
          ("epsilons", po::value<std::string>()->default_value(""), "epsilon values for simplification");

  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cerr << desc << std::endl;
    return 1;
  }

  if (vm.count("dataset") == 0 || vm.count("queries") == 0 ||
      vm.count("algorithm") == 0) {
    std::cerr << desc << std::endl;
    return 1;
  }

  return 0;
}

extern std::vector<std::vector<size_t>> tl_results;
#pragma omp threadprivate(tl_results)
std::vector<std::vector<size_t>> tl_results;
extern std::vector<std::pair<double, id_pair_t>> tl_scored_results;
#pragma omp threadprivate(tl_scored_results)
std::vector<std::pair<double, id_pair_t>> tl_scored_results;

std::vector<double> get_epsilon_vec(std::string &epsilon_str) {
  std::vector<double> epsilons;
  if (epsilon_str.size() > 0) {
    std::vector<std::string> tokens;
    split(tokens, epsilon_str, boost::algorithm::is_any_of(","));
    for(const auto & token : tokens) {
      epsilons.push_back(atof(token.c_str()));
    }
  }
  return epsilons;
}

inline bool find_vec(const std::vector<size_t> & v, size_t x) {
  for (auto e : v) {
    if (x == e) return true;
  }
  return false;
}

/**
 * Estimate the evaluation thresholds. 
 */
template<typename Point>
double
estimate_thresholds(const std::vector<Point> & queries,
                    const std::vector<Point> & dataset,
                    const double frac_pairs_to_eval,
                    const size_t sample_size,
                    std::vector<std::pair<double, std::pair<size_t, size_t>>> & scored_pairs,
                    Xorshift1024star & rnd) {
  const size_t num_buckets = 1000;
  std::vector<size_t> hist(num_buckets);

  // Sample the histogram
  size_t nonzero = 0;
  for (const auto & p : scored_pairs) {
    const double score = p.first;
    if (score > 0.0) {
      nonzero++;
      size_t bucket = (size_t) (score * num_buckets);
      if (bucket == num_buckets) bucket--;
      hist[bucket]++;
    }
  }

  const size_t count_threshold = frac_pairs_to_eval * nonzero;
  /* printf("Count threshold is %lu\n", count_threshold); */
  size_t i = 0;
  size_t count = 0;
  for (; count<count_threshold; i++) {
    count += hist[i];
    /* printf("Count %lu\n", count); */
  }
  const double upper_score = i / ((double) num_buckets);

  std::cout << "Upper score for evaluation " << upper_score << std::endl;

  return upper_score;
}

template<typename Point>
void idxs_to_ids(const std::vector<Point> & queries,
                 const std::vector<Point> & dataset,
                 std::vector<std::pair<double, std::pair<size_t, size_t>>> & sr) {
  for (auto & p : sr) {
    size_t query_idx = p.second.first;
    size_t data_idx = p.second.second;
    size_t query_id = queries[query_idx].id;
    size_t data_id = dataset[data_idx].id;
    p.second.first = query_id;
    p.second.second = data_id;
  }
}

template <typename Point>
int main_template(const boost::program_options::variables_map &vm) {
  START_TIMER(total);
  printf("Using %d threads\n", omp_get_max_threads());

  std::string dataset_path = vm["dataset"].as<std::string>();
  std::string queries_path = vm["queries"].as<std::string>();
  std::string algorithm = vm["algorithm"].as<std::string>();
  std::string baseline_basedir = vm["baseline"].as<std::string>();
  std::string eval_dist_str =  vm["eval-dist"].as<std::string>();
  bool eval_dist =  eval_dist_str.compare("yes") == 0;
  double verify_fraction = 0;
  if (eval_dist) {
    printf("Evaluating distance\n");
  } else {
    try {
      verify_fraction = std::stof(eval_dist_str);
    } catch (...) {
      printf("Could not parse %s into a double, using 0.0 instead\n", eval_dist_str.c_str()); 
    }
    printf("Verifying the %f lowest-proability collisions\n", verify_fraction);
  }
  double range = vm["range"].as<double>();
  EXPERIMENT_TAG("dataset", dataset_path);
  EXPERIMENT_TAG("queries", queries_path);
  EXPERIMENT_TAG("algorithm", algorithm);
  EXPERIMENT_TAG("range", range);
  EXPERIMENT_TAG("eval-dist", vm["eval-dist"].as<std::string>());

  uint64_t seed;
  if (vm.count("seed")) {
    seed = vm["seed"].as<uint64_t>();
  } else {
    std::random_device rd;
    seed = rd();
  }
  printf("Using random seed %lu\n", seed);
  EXPERIMENT_TAG("seed", std::to_string(seed));
  Splitmix64 seeder(seed);
  Xorshift1024star rnd(seeder.next());

  std::string epsilon_str = vm["epsilons"].as<std::string>();
  EXPERIMENT_TAG("epsilons", epsilon_str);
  std::vector<double> epsilons = get_epsilon_vec(epsilon_str);

  START_TIMER(data_loading);
  std::vector<Curve<Point>> dataset = load_dataset_msgpack<Point>(dataset_path);
  std::vector<Curve<Point>> queries = load_dataset_msgpack<Point>(queries_path);
  STOP_TIMER_V(data_loading);
  const size_t n_queries = queries.size();
#define RESULT_IDX(query_idx, data_idx) query_idx*n_queries + data_idx
  std::vector<std::pair<double, std::pair<size_t, size_t>>> scored_results;

  // to be used with verification
  std::function<void(size_t, size_t)> positive_callback = [&](size_t query_idx, size_t data_idx) {
    tl_results[query_idx].push_back(data_idx);
    return;
  };

  // to be used with partial verification or no verification
  std::function<void(size_t, size_t, double)> partial_verification_callback =
    [&](size_t query_idx, size_t data_idx, double score) {
      // reports the IDs, not the indices
      size_t query_id = queries[query_idx].id;
      size_t data_id = dataset[data_idx].id;
      scored_results[RESULT_IDX(query_idx, data_idx)] = 
        std::make_pair(score, std::make_pair(query_idx, data_idx));
    };

  if (algorithm.compare("exact") == 0) {
#pragma omp parallel
    {
      tl_results.resize(n_queries);
    };
    query_continuous_exact(dataset, queries, range, epsilons, rnd, positive_callback);
  } else if (algorithm.compare("basic") == 0) {
#pragma omp parallel
    {
      tl_results.resize(n_queries);
    };
    query_continuous_basic(dataset, queries, range, rnd, positive_callback);
  } else if (algorithm.compare("hash") == 0) {
    scored_results.resize(queries.size() * dataset.size());
    double resolution_factor = 4*Point::dimensions();
    if (vm.count("resolution_factor") > 0){
      resolution_factor = vm["resolution-factor"].as<double>();
    }
    std::cout << "Resolution factor " << resolution_factor << std::endl;
    size_t k = vm["k"].as<size_t>();
    size_t L = vm["L"].as<size_t>();
    if (L==0) {
      // Set L based on k as 0.5^{-k}
      L = 1 << k;
      printf("L set to %lu\n", L);
    }
    EXPERIMENT_TAG("k", k);
    EXPERIMENT_TAG("L", L);
    EXPERIMENT_TAG("eval_dist", eval_dist);
    EXPERIMENT_TAG("resolution_factor", resolution_factor);
    if (eval_dist) {
      query_hash(dataset, queries, range, epsilons, k, L, 
                 resolution_factor, rnd, positive_callback);
    } else {
      query_hash_no_eval(dataset, queries, range, epsilons, k, L, 
                         resolution_factor, rnd, partial_verification_callback);
      if (verify_fraction > 0) {
        START_TIMER(threshold_estimation);
        const double upper_bound =
          estimate_thresholds(queries, dataset, verify_fraction, 10000, scored_results, rnd);
        STOP_TIMER_V(threshold_estimation);
        START_TIMER(distance_verification);
        size_t verified_cnt = 0;
#pragma omp parallel for
        for (size_t i=0; i<scored_results.size(); i++) {
          auto & p = scored_results[i];
          const double score = p.first;
          if (score > 0.0 && score <= upper_bound) {
#pragma omp atomic
            verified_cnt++;
            // Verify this pair because is in the undecided zone
            const size_t query_idx = p.second.first;
            const size_t data_idx = p.second.second;
            if (continuous_frechet_distance_predicate(
                  queries[query_idx], dataset[data_idx], range, epsilons)) {
              p.first = 1.0;
            } else {
              p.first = 0.0;
            }
          }
        }
        STOP_TIMER_V(distance_verification);
        printf("Verified pairs %lu\n", verified_cnt);
        EXPERIMENT_APPEND("auto_estimation", {{"upper_bound", upper_bound},
                                              {"num_verified", verified_cnt}});
      }
    }
  } else {
    throw std::logic_error("Unknown algorithm");
  }

  STOP_TIMER_V(total);
  report_counters();

  printf("Adding results to the experiment\n");
  if (eval_dist) {
    std::vector<std::vector<size_t>> result;
    result.resize(n_queries);
#pragma omp parallel
    {
#pragma omp critical
      {
        for (size_t i=0; i<n_queries; i++) {
          result[i].insert(result[i].end(), tl_results[i].begin(), tl_results[i].end());
        }
      }
    }

    auto pairs = result_as_pairs(result, queries, dataset);
    if (algorithm.compare("exact") == 0 || algorithm.compare("basic") == 0) {
      // then dump the baseline
      printf("Dumping the baseline\n");
      dump_baseline(pairs, baseline_basedir, queries_path, dataset_path, range);
    }

    printf("Reading baseline\n");
    auto baseline = load_baseline(baseline_basedir, queries_path, dataset_path, range);

    double recall_score = recall(pairs, baseline.pairs);
    double precision_score = precision(pairs, baseline.pairs);
    print_false_positives(pairs, baseline.pairs, queries, dataset);
    double auc_score = auc(0.0, recall_score);
    size_t output_size = pairs.size();
    double total_time = GET_WALL_SUM(std::string("total"));
    double ms_per_pair = total_time / output_size;
    printf("Milliseconds per output pair: %f\n", ms_per_pair);
    EXPERIMENT_APPEND("scores", {{"auc", auc_score}, {"precision", precision_score},{"ms_per_pair", ms_per_pair}});
    printf("Recall of %f\n", recall_score);
    printf("Precision of %f\n", precision_score);
    printf("AUC of %f\n", auc_score);
  } else {
    idxs_to_ids(queries,dataset, scored_results);
    printf("Reading baseline\n");
    auto baseline = load_baseline(baseline_basedir, queries_path, 
                                  dataset_path, range);
    size_t total_pairs = dataset.size() * queries.size();
    size_t output_size = scored_results.size();
    double total_time = GET_WALL_SUM(std::string("total"));
    double ms_per_pair = total_time / output_size;
    /* auto roc_curve = roc(scored_results, baseline.pairs, total_pairs); */
    /* auto pr_curve = precision_recall_curve(scored_results, baseline.pairs, total_pairs); */
    double all_recall = recall(scored_results, baseline.pairs);
    double all_precision = precision(scored_results, baseline.pairs);
    double all_fscore = 2*all_precision*all_recall / (all_precision + all_recall);
    // TODO uncomment the following line
    /* report_scores(scored_results, baseline.pairs); */
    if (vm.count("dump-solution") > 0) {
      std::string output_path = vm["dump-solution"].as<std::string>();
      std::cout << "Dumping the solution to " << output_path << std::endl;
      dump_result(scored_results, baseline.pairs, output_path);
    }
    /* for(auto e : roc_curve) { */
    /*   EXPERIMENT_APPEND("roc-curve", {{"fpr", e.first}, */
    /*                                   {"tpr", e.second}}); */
    /* } */
    /* printf("  -- Recall Precision\n"); */
    /* for(auto e : pr_curve) { */
    /*   /1* printf("  -- %f %f\n", e.first, e.second); *1/ */
    /*   EXPERIMENT_APPEND("pr-curve", {{"precision", e.second}, */
    /*                                  {"recall", e.first}}); */
    /* } */
    /* double auc_score = auc(roc_curve); */
    /* double auc_pr_score = auc(pr_curve); */
    /* printf("AUC ROC score %f\n", auc_score); */
    /* printf("AUC PR score %f\n", auc_pr_score); */
    printf("precision score %f\n", all_precision);
    printf("recall score %f\n", all_recall);
    printf("F-score %f\n", all_fscore);
    EXPERIMENT_APPEND("scores", {//{"auc_pr_roc", auc_pr_score},
                                 {"fscore_all", all_fscore},
                                 {"precision_all", all_precision}, 
                                 {"recall_all", all_recall},
                                 {"ms_per_pair", ms_per_pair}});
  }

  printf("Adding timers to the experiment\n");
  EXPERIMENT_RECORD_PROFILE();
  printf("Saving experiment\n");
  START_TIMER(save_experiment);
  EXPERIMENT_SAVE();
  STOP_TIMER_V(save_experiment);

  return 0;
}


int main(int argc, char **argv) {
#ifndef NDEBUG
  printf("WARNING: The code is running in debug mode (or at least with assertions enabled)\n");
#endif

  boost::program_options::variables_map vm;
  int ret = parse_opts(argc, argv, vm);
  if (ret != 0) {
    return ret;
  }
  std::string dataset_path = vm["dataset"].as<std::string>();
  std::string queries_path = vm["queries"].as<std::string>();
  size_t dimensions = read_dataset_dimensions(dataset_path);
  std::cout << "Dataset of points of dimensionality " 
            << dimensions
            << " queries with points of dimensionality " 
            << read_dataset_dimensions(dataset_path)
            << std::endl;

  printf("Running version %s\n", g_GIT_SHA1);

  switch(dimensions) {
    case 1: 
      return main_template<Point1D>(vm);
    case 2: 
      return main_template<Point2D>(vm);
    default:
      std::cerr << "Unsupported number of dimensions " << dimensions << std::endl;
      return -1;
  }
}

