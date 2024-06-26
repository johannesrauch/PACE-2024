cmake_minimum_required(VERSION 3.15)

project(pace LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Release)
if(PACE_DEBUG_PRINT)
    add_compile_definitions(PACE_DEBUG_PRINT)
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)
file(GLOB PACE_SOURCES ${CMAKE_SOURCE_DIR}/src/*/*.cpp)

# Uncomment the following line if you do not have highs installed
add_subdirectory(highs)
find_package(HIGHS REQUIRED)
find_package(Threads REQUIRED)
add_subdirectory(src)

add_executable(weberknecht ${CMAKE_SOURCE_DIR}/src/weberknecht.cpp)
target_link_libraries(weberknecht oscm)
add_executable(weberknecht_h ${CMAKE_SOURCE_DIR}/src/weberknecht_h.cpp)
target_link_libraries(weberknecht_h oscm)

if(PACE_BUILD_TESTS)
    add_subdirectory(test)
endif()
