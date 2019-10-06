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

#ifndef EXPERIMENT_REPORTER_H
#define EXPERIMENT_REPORTER_H

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/variant.hpp>
#include <boost/variant/get.hpp>
#include <exception>
#include <fstream>
#include <map>
#include <set>
#include "msgpack.hpp"
#include "types.h"

typedef boost::variant<size_t, int, double, std::string, std::vector<size_t>, bool> element_value_t;

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
    namespace adaptor {

    template <>
    struct pack<element_value_t> {
        template <typename Stream>
        msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, element_value_t const& v) const {
            switch (v.which()) {
                case 0:  // size_t
                {
                    o.pack_unsigned_long_long(boost::get<size_t>(v));
                    break;
                }
                case 1:  // int
                {
                    o.pack_long_long(boost::get<int>(v));
                    break;
                }
                case 2:  // double
                {
                    o.pack_double(boost::get<double>(v));
                    break;
                }
                case 3:  // std::string
                {
                    o.pack(boost::get<std::string>(v));
                    break;
                }
                case 4:  // std::vector<size_t>
                {
                    o.pack(boost::get<std::vector<size_t>>(v));
                    break;
                }
                case 5:  // bool
                {
                    if (boost::get<bool>(v))
                        o.pack_true();
                    else
                        o.pack_false();
                    break;
                }
            }
            return o;
        }
    };

    template <>
    struct pack<boost::posix_time::ptime> {
        template <typename Stream>
        msgpack::packer<Stream>& operator()(msgpack::packer<Stream>& o, boost::posix_time::ptime const& v) const {
            std::string s = boost::posix_time::to_iso_extended_string(v);
            o.pack(s);
            return o;
        }
    };

    }  // namespace adaptor
}  // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
}  // namespace msgpack

class ExperimentReporter {
    typedef std::map<std::string, element_value_t> table_row_t;
    typedef std::vector<table_row_t> table_t;

  private:
    ExperimentReporter() : date(boost::posix_time::second_clock::local_time()) {}

  public:
    static ExperimentReporter& get_instance() {
        static ExperimentReporter instance;
        return instance;
    }

    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    void set_results(std::vector<std::vector<std::pair<size_t, double>>>* res) {
        results = res;
    }
    // void set_results(std::vector<std::vector<size_t>>* res) {
    //     results = res;
    // }

    void tag(const std::string& key, const element_value_t& val) {
        tags[key] = val;
    }

    void append(const std::string& table, table_row_t row) {
        if (tables.find(table) == tables.end()) {
            tables[table] = table_t();
        } else {
            table_row_t first_row = tables[table].front();
            std::set<std::string> keys;
            for (auto it = first_row.begin(); it != first_row.end(); ++it) {
                keys.insert(it->first);
            }
            if (row.size() != keys.size()) {
                throw std::runtime_error("New row and existing ones have a different number of columns");
            }
            for (auto it = row.begin(); it != row.end(); ++it) {
                if (keys.find(it->first) == keys.end()) {
                    throw std::runtime_error("New row has different headings from the existing ones");
                }
            }
        }
        tables[table].push_back(row);
    }

    void save_msgpack(std::ostream& out) {
        msgpack::packer<std::ostream> pk(&out);
        pk.pack_map(3);  // The three top level fields

        // First field of map
        pk.pack(std::string("date"));
        pk.pack(boost::posix_time::to_iso_extended_string(date));

        // Second field of map
        pk.pack(std::string("tags"));
        pk.pack(tags);

        // Third field in the map, the tables. We do each table by streaming
        pk.pack(std::string("tables"));
        pk.pack_map(tables.size() + 1);  // the +1 is to store the table of the QueryResults
        for (auto tables_it = tables.begin(); tables_it != tables.end(); ++tables_it) {
            std::string tname = tables_it->first;
            pk.pack(tname);
            table_t table = tables_it->second;
            // Stream the packing of rows
            pk.pack_array(table.size());
            for (auto row_it = table.begin(); row_it != table.end(); ++row_it) {
                pk.pack(*row_it);
            }
        }

        /// FIXME: Report the results
        pk.pack(std::string("results"));
        if (results == nullptr) {
            pk.pack_array(0);
        } else {
            const size_t n_results = results->size();
            pk.pack_array(n_results);
            for (size_t i = 0; i < n_results; i++) {
                pk.pack_map(2);
                pk.pack(std::string("query"));
                pk.pack(i);
                pk.pack(std::string("matches"));
                pk.pack_array(results->at(i).size());
                for (auto x : results->at(i)) {
                    pk.pack(x);
                }
            }
        }
    }

    void save_msgpack(const std::string& path) {
#ifdef EXPERIMENT_COMPRESS
        std::ofstream out_stream(path, std::ios_base::app | std::ios_base::out);
        boost::iostreams::filtering_ostream out;
        out.push(boost::iostreams::bzip2_compressor());
        out.push(out_stream);
        save_msgpack(out);
#else
        std::ofstream out(path, std::ios_base::app | std::ios_base::out);
        save_msgpack(out);
        out.close();
#endif
    }

    void save_msgpack() {
        // std::string path = boost::posix_time::to_iso_string(date);
#ifdef EXPERIMENT_COMPRESS
        std::string path = "results.mpk.bz2";
#else
        std::string path = "results.mpk";
#endif
        save_msgpack(path);
    }

  private:
    boost::posix_time::ptime date;
    std::map<std::string, element_value_t> tags;
    std::map<std::string, table_t> tables;
    // Modified by Shunsuke Kanda (https://github.com/kampersanda)
    std::vector<std::vector<std::pair<size_t, double>>>* results = nullptr;
    // std::vector<std::vector<size_t>>* results = nullptr;
};

#define EXPERIMENT_SET_RESULTS(results_ptr) ExperimentReporter::get_instance().set_results(results_ptr)

#define EXPERIMENT_APPEND(table, ...) ExperimentReporter::get_instance().append(table, __VA_ARGS__)

#define EXPERIMENT_TAG(key, val) ExperimentReporter::get_instance().tag(std::string(key), val)

#define EXPERIMENT_SAVE() ExperimentReporter::get_instance().save_msgpack()

// Added by Shunsuke Kanda (https://github.com/kampersanda)
#define EXPERIMENT_SAVE_WITH_PATH(path) ExperimentReporter::get_instance().save_msgpack(std::string(path))

#endif /* EXPERIMENT_REPORTER_H */
