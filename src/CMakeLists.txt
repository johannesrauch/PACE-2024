add_library(oscm STATIC ${PACE_SOURCES})
target_link_libraries(oscm highs::highs)