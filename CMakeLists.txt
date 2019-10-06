cmake_minimum_required(VERSION 3.0)
project(frechet_simsearch)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# C++17 compiler check
if(CMAKE_CXX_COMPILER MATCHES ".*clang.*" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_COMPILER_IS_INTEL 1)
endif()
if((CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 7.0) OR (CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.0))
  message(FATAL_ERROR "Your C++ compiler does not support C++17. Please install g++ 7.0 (or greater) or clang 4.0 (or greater)")
else()
  message(STATUS "Compiler is recent enough to support C++17.")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-fdiagnostics-color=always)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_compile_options(-fcolor-diagnostics)
endif()

set(BUILTIN_POPCNT 0)

if(DISABLE_SSE4_2)
  message(STATUS "sse4.2 disabled")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # Use /proc/cpuinfo to get the information
  file(STRINGS "/proc/cpuinfo" _cpuinfo)
  if(_cpuinfo MATCHES "(sse4_2)|(sse4a)")
    set(BUILTIN_POPCNT 1)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  execute_process(COMMAND sysctl -n machdep.cpu.features OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_cpuinfo MATCHES "SSE4.2")
    set(BUILTIN_POPCNT 1)
  endif()
endif()

if(BUILTIN_POPCNT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
endif()

# Dependencies
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED
  COMPONENTS system filesystem iostreams timer chrono program_options date_time)
include_directories(${Boost_INCLUDE_DIRS})

# SDSL (Please modify these depending on your env.)
link_directories($ENV{HOME}/lib)
include_directories($ENV{HOME}/include)

find_package(BZip2)
if(BZIP2_FOUND)
  add_definitions(-DEXPERIMENT_COMPRESS=1)
else()
  message("-- BZip2 not found, experimental results will not be compressed")
endif()

set(GCC_WARNINGS "-Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -pthread ${GCC_WARNINGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -march=native -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -O0 -g -DDEBUG")

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

message(STATUS "BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX_FLAGS are ${CMAKE_CXX_FLAGS}")
message(STATUS "CXX_FLAGS_DEBUG are ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CXX_FLAGS_RELEASE are ${CMAKE_CXX_FLAGS_RELEASE}")

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    include_directories(/usr/local/include /usr/local/opt/openssl/include)
    link_directories(/usr/local/lib /usr/local/opt/openssl/lib)
endif()

add_subdirectory(external)
include_directories(external)
include_directories(include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(EXE_SOURCE_FILES
    analyze_distance.cpp
    analyze_array_perf.cpp
    analyze_array_score.cpp
    analyze_fresh_perf.cpp
    analyze_fresh_score.cpp
    analyze_hmsearch_perf.cpp
    analyze_mtrie_perf.cpp
    make_range_groundtruth.cpp
    run_exact_search.cpp
)

foreach(EXE_SOURCE_FILE ${EXE_SOURCE_FILES})
    MESSAGE("TARGET:" ${EXE_SOURCE_FILE})
    get_filename_component(EXE_PREFIX ${EXE_SOURCE_FILE} NAME_WE)
    add_executable(${EXE_PREFIX} src/${EXE_SOURCE_FILE})
    target_link_libraries(${EXE_PREFIX} curvedist sdsl ${Boost_LIBRARIES} ${BZIP2_LIBRARIES})
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        target_link_libraries(${EXE_PREFIX} omp)
    endif()
endforeach()
