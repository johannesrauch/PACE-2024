
set(PACE_ALL_TESTS
bipartite_graph
branch_and_cut
crossings_utils
exact_solver
heuristics
highs_mip
index
input
matrix
parse_input
topological_sort
transitive_hull
)

foreach(PACE_TEST IN LISTS PACE_ALL_TESTS)
    message(STATUS "add_executable(${PACE_TEST} ${PACE_TEST}.cpp)")
    add_executable(${PACE_TEST} ${PACE_TEST}.cpp)
    target_link_libraries(${PACE_TEST} highs::highs)
endforeach()

enable_testing()

foreach(PACE_TEST IN LISTS PACE_ALL_TESTS)
    add_test(
        NAME ${PACE_TEST} 
        COMMAND $<TARGET_FILE:${PACE_TEST}> 
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test
    )
endforeach()
