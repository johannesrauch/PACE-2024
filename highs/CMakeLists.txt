project(LOAD_HIGHS LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    highs
    GIT_REPOSITORY "https://github.com/ERGO-Code/HiGHS.git"
    GIT_TAG        "v1.7.0"
)
set(FAST_BUILD ON CACHE INTERNAL "Fast Build")

FetchContent_MakeAvailable(highs)
