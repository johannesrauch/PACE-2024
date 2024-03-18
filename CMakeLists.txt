# set minimum cmake version
cmake_minimum_required(VERSION 3.15)

# project name and language
project(pace2024 LANGUAGES CXX)

# require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LINK_FLAGS -lhighs -lstdc++ -lm)

add_subdirectory(test)