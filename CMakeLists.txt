# set minimum cmake version
cmake_minimum_required(VERSION 3.15)

# project name and language
project(pace LANGUAGES CXX)

# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# include and link directories
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/highs/include/highs)
link_directories(${CMAKE_SOURCE_DIR}/highs/lib)
#find_library(highs REQUIRED)

# build
add_subdirectory(src)
add_subdirectory(test)
