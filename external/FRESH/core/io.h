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

#ifndef CURVEDIST_IO_H
#define CURVEDIST_IO_H

#include "types.h"
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "timer.h"
#include <stdio.h>

std::map<std::string, std::string> get_file_metadata(const std::string & path);

size_t read_dataset_dimensions(std::string & path);

std::vector< Curve<Point2D> > load_dataset_giscup(const std::string & path);

template <typename Point>
std::vector< Curve<Point> > load_dataset_msgpack(const std::string & path) {
  bool has_metadata = get_file_metadata(path).count("dimensions") > 0;
  std::vector< Curve<Point> > dataset;
  std::ifstream input_stream(path.c_str(), std::ifstream::in);
  if (!input_stream.good()) {
    throw std::logic_error("Could not open file");
  }

  boost::iostreams::filtering_istream fis;
  if (boost::algorithm::ends_with(path, ".bz2")) {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(input_stream);

  msgpack::unpacker pac; // Object to handle the deserialization on the go.

  size_t try_read_size = 64*1024;
  
  size_t cnt = 0;
  while (!fis.eof()) {
    pac.reserve_buffer(try_read_size);
    fis.read(pac.buffer(), try_read_size);
    size_t actual_read_bytes = fis.gcount();
    pac.buffer_consumed(actual_read_bytes);

    msgpack::object_handle result;
    while(pac.next(result)) {
      if (cnt > 0 || !has_metadata) {
        Curve<Point> c;
        result.get().convert(c);
        c.fix_prefix_lengths();
        dataset.push_back(c);
      } else {
        std::cout << "Skipping metadata " << result.get() << std::endl;
      }
      cnt++;
    }
  }  
  return dataset;
}

std::vector< Curve<Point2D> > load_dataset(const std::string & path);


#endif //CURVEDIST_IO_H
