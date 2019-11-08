#pragma once

#include <cxxabi.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <regex>
#include <string>
#include <vector>

#include <boost/algorithm/string/classification.hpp>  // is_any_of
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>

#ifdef __SSE4_2__
#include <xmmintrin.h>
#endif

#include <tinyformat.h>

namespace frechet_simsearch {

static const uint32_t DOMAIN_BITS = 8;
static const char* NO_COLOR = "\033[0;0m";

template <typename... Args>
inline void errorfln(const char* fmt, const Args&... args) {
    const char* RED_COLOR = "\033[0;31m";
    std::cerr << RED_COLOR << "ERROR: ";
    tfm::format(std::cerr, fmt, args...);
    std::cerr << NO_COLOR << std::endl;
}
template <typename... Args>
inline void warnfln(const char* fmt, const Args&... args) {
    const char* YELLOW_COLOR = "\033[0;33m";
    std::cerr << YELLOW_COLOR << "WARNING: ";
    tfm::format(std::cerr, fmt, args...);
    std::cerr << NO_COLOR << std::endl;
}
template <typename... Args>
inline void reportfln(const char* fmt, const Args&... args) {
    const char* GREEN_COLOR = "\033[0;32m";
    std::cout << GREEN_COLOR;
    tfm::format(std::cout, fmt, args...);
    std::cout << NO_COLOR << std::endl;
}

struct Configure {
    uint32_t concatenations = 0;
    uint32_t repeatations = 0;
    uint32_t buckets = 0;  // # of buckets (of Multi-Index)
    uint32_t reduction_threshold = 0;  // threshould of node reduction
    double resolution = 0.0;  // of grid
    double resolution_factor = 0.0;  // of grid
    uint64_t seed = 0;

    std::string to_str() const {
        return tfm::format("K_%d-L_%d-B_%d-T_%d-R_%g-F_%g-S_%d",  //
                           concatenations, repeatations, buckets, reduction_threshold, resolution, resolution_factor,
                           seed);
    }
    void show(std::ostream& os) const {
        tfm::format(os, "[configure]\n");
        tfm::format(os, "- concatenations: %d\n", concatenations);
        tfm::format(os, "- repeatations: %d\n", repeatations);
        tfm::format(os, "- buckets: %d\n", buckets);
        tfm::format(os, "- reduction_threshold: %d\n", reduction_threshold);
        tfm::format(os, "- resolution: %g\n", resolution);
        tfm::format(os, "- resolution_factor: %g\n", resolution_factor);
        tfm::format(os, "- seed: %d\n", seed);
    }
};

inline std::string normalize_filepath(std::string filepath) {
    std::replace(filepath.begin(), filepath.end(), '/', '_');
    std::replace(filepath.begin(), filepath.end(), '.', '_');
    std::replace(filepath.begin(), filepath.end(), ':', '_');
    return filepath;
}

inline std::string getnow() {
    char data[64];
    std::time_t t = std::time(NULL);
    std::strftime(data, sizeof(data), "%Y_%m_%d_%H_%M_%S", std::localtime(&t));
    return std::string(data);
}

inline bool exists_path(const std::string path) {
    boost::filesystem::path p(path);
    return boost::filesystem::exists(p);
}

inline void make_directory(const std::string dir) {
    boost::filesystem::path path(dir);
    if (!boost::filesystem::exists(path)) {
        if (!boost::filesystem::create_directories(path)) {
            errorfln("unable to create output directory %s", dir);
            exit(1);
        }
    }
}

inline std::tuple<uint32_t, uint32_t, uint32_t> parse_range(std::string range_str) {
    std::vector<std::string> result;
    boost::algorithm::split(result, range_str, boost::is_any_of(":"));

    if (result.size() == 1) {
        uint32_t max = std::stoi(result[0]);
        return {0, max, 1};
    }
    if (result.size() == 2) {
        uint32_t min = std::stoi(result[0]);
        uint32_t max = std::stoi(result[1]);
        return {min, max, 1};
    }
    if (result.size() == 3) {
        uint32_t min = std::stoi(result[0]);
        uint32_t max = std::stoi(result[1]);
        uint32_t stp = std::stoi(result[2]);
        return {min, max, stp};
    }

    errorfln("Invalid format of range string %s", range_str);
    exit(1);
}

class ProgressPrinter {
  public:
    ProgressPrinter(size_t end) {
        const size_t step = end / 10;
        for (size_t i = 1; i <= 9; ++i) {
            m_stations[i - 1] = i * step;
        }
        m_stations[9] = end;
        m_stations[10] = size_t(-1);
        std::cout << "# " << std::flush;
    }

    void operator()(size_t i) {
        while (m_stations[m_cursor] <= i) {
            std::cout << "* " << std::flush;
            ++m_cursor;
        }
        if (m_cursor == 10) {
            std::cout << "done!!" << std::endl;
        }
    }

  private:
    size_t m_stations[11];
    size_t m_cursor = 0;
};

}  // namespace frechet_simsearch
