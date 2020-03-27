cmake_minimum_required(VERSION 3.4)
project(Utils)
set(CMAKE_BUILD_TYPE "Release")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(.)
set(CMAKE_CXX_STANDARD 11)

add_library(utils SHARED Utils.cpp)
add_subdirectory(test)

install(TARGETS utils LIBRARY DESTINATION lib)
install(FILES Utils.h DESTINATION include/lyxutils/)
install(FILES UtilsConfig.cmake DESTINATION lib/cmake/Utils/)