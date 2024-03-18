include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(/usr/local/include/highs)

set(ALL_TESTS
bipartite_graph
branch_and_cut_tiny
branch_and_cut_medium
crossings
debug_printf
heuristics
input
lp_wrapper
matrix
median_heuristic
shift_heuristic
topological_sort
)

foreach(T IN LISTS ALL_TESTS)
    message(STATUS "add_executable(${T} ${T}.cpp)")
    add_executable(${T} ${T}.cpp)
    target_link_libraries(${T} highs)
endforeach()

enable_testing()

foreach(T IN LISTS ALL_TESTS)
    add_test(
        NAME ${T} 
        COMMAND $<TARGET_FILE:${T}> 
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    )
endforeach()
