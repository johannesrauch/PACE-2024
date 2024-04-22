
set(ALL_TESTS
bipartite_graph
branch_and_cut
crossings_utils
heuristics
input
matrix
parse_input
topological_sort
transitive_hull
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
