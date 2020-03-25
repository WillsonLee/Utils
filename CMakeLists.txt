cmake_minimum_required(VERSION 3.15)
project(Utils)
set(CMAKE_BUILD_TYPE "Release")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(.)
set(CMAKE_CXX_STANDARD 11)

add_library(utils SHARED Utils.cpp)

INSTALL(TARGETS utils LIBRARY DESTINATION /usr/local/lib/)
INSTALL(FILES Utils.h DESTINATION /usr/local/include/lyxutils/)

add_subdirectory(test)