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

#include "io.h"

Curve<Point2D> load_curve_giscup(const std::string & path) {
  std::ifstream input(path);
  if (!input.good()) {
    throw std::runtime_error("File not found: " + path);
  }

  std::vector<Point2D> curve;
  size_t id = std::numeric_limits<size_t>::max();
  bool parse_id = true;
  
  double last_x = std::numeric_limits<double>::infinity();
  double last_y = std::numeric_limits<double>::infinity();

  std::string line;
  std::getline(input, line); // Skip first line, which contains metadata
  while(std::getline(input, line)) {
    std::stringstream sstream(line);
    double x, y;
    sstream >> x >> y;
    if (parse_id) {
      size_t i;
      sstream >> i >> id; 
    }
    if (last_x != x || last_y != y) {
      last_x = x;
      last_y = y;
      curve.emplace_back(x, y);
    }
  }
  return Curve<Point2D>(id, curve);
}

std::vector< Curve<Point2D> > load_dataset_giscup(const std::string & path) {
  std::ifstream input(path);
  if (!input.is_open()) {
    throw std::logic_error("Error opening " + path);
  }

  std::vector< Curve<Point2D> > curves;
  
  std::string line;
  while (std::getline(input, line)) {
    Curve<Point2D> c = load_curve_giscup(line);
    if (c.points.size() > 0) {
      curves.emplace_back(c);
    }
  }

  return curves;
}

std::vector< Curve<Point2D> > load_dataset(const std::string & path){
  namespace ba = boost::algorithm;
  if (ba::ends_with(path, ".txt")) {
    std::cerr << "Loading dataset in GisCup format" << std::endl;
    return load_dataset_giscup(path);
  } else {
    std::cerr << "Loading dataset in MsgPack format" << std::endl;
    return load_dataset_msgpack<Point2D>(path);
  }
}

std::map<std::string, std::string>
get_file_metadata(const std::string & path) {
  const std::map<std::string, std::string> empty;
  namespace ba = boost::algorithm;
  if (ba::ends_with(path, ".txt")) {
    return empty;
  } else {
    std::ifstream input_stream(path.c_str(), std::ifstream::in);
    if (!input_stream.good()) {
      throw std::logic_error("Could not open file");
    }

    boost::iostreams::filtering_istream fis;
    if (boost::algorithm::ends_with(path, ".bz2")) {
      fis.push(boost::iostreams::bzip2_decompressor());
    }
    fis.push(input_stream);

    msgpack::unpacker pac;

    size_t try_read_size = 1024;
    
    while (!fis.eof()) {
      pac.reserve_buffer(try_read_size);
      fis.read(pac.buffer(), try_read_size);
      size_t actual_read_bytes = fis.gcount();
      pac.buffer_consumed(actual_read_bytes);

      msgpack::object_handle result;
      while(pac.next(result)) {
        try {
          std::map<std::string, std::string> metadata;
          result.get().convert(metadata);
          return metadata;
        } catch(msgpack::type_error) {
          // The map we just read represents a point, so there is no metadata
          // Since a point is a heterogeneous map, then we have to catch 
          // an exception
          return empty;
        }
      }
    }
    return empty; // This point is actually never reached
  }
}

size_t read_dataset_dimensions(std::string & path) {
  auto meta = get_file_metadata(path);
  if (meta.count("dimensions") > 0) {
    return std::stoi(meta["dimensions"]);
  } else {
    return 2; // 2 is the number of dimensions of the legacy datasets without metadata
  }
}


