## Core build

# Commented out by Shunsuke Kanda (https://github.com/kampersanda)
# configure_file(
#   "${CMAKE_CURRENT_SOURCE_DIR}/git_info.cpp.in"
#   "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp" @ONLY)

set(CURVEDIST_LIB_HEADER_FILES
  prelude.h
  # Commented out by Shunsuke Kanda (https://github.com/kampersanda)
  # git_info.h
  frechet.h
  types.h
  io.h
  timer.h
  threading.h
  hash.h
  query.h
  stats.h
  densify.h
  logging.h
  experiment_reporter.h
  babr/geometry_basics.hpp
  babr/frechet_distance.hpp
  babr/frechet_distance2.hpp)

set(CURVEDIST_LIB_SOURCE_FILES
  types.cpp
  frechet.cpp
  io.cpp
  densify.cpp
  hash.cpp
  timer.cpp
  query.cpp
  stats.cpp
  babr/frechet_distance.cpp
  babr/geometry_basics.cpp
  # Commented out by Shunsuke Kanda (https://github.com/kampersanda)
  # "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp"
  scores.h)

add_library(curvedist STATIC
        ${CURVEDIST_LIB_HEADER_FILES}
        ${CURVEDIST_LIB_SOURCE_FILES})

# Commented out by Shunsuke Kanda (https://github.com/kampersanda)
# add_executable(cdist main.cpp ${CURVEDIST_LIB_HEADER_FILES})
# target_link_libraries(cdist curvedist
#   ${Boost_LIBRARIES}
#   ${BZIP2_LIBRARIES})
