cmake_minimum_required(VERSION 3.15)
project(pace LANGUAGES CXX)

# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Release)

# include and link directories
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/highs/include/highs)

set(HIGHS_DIR ${CMAKE_SOURCE_DIR}/lib/cmake/highs)
find_package(HIGHS REQUIRED)
find_package(Threads REQUIRED)

add_executable(weberknecht src/weberknecht.cpp)
target_link_libraries(weberknecht highs::highs)

if(PACE_BUILD_TESTS)
    add_subdirectory(test)
endif()
