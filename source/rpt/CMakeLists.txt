file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/source/rpt/*.cc" "${CMAKE_SOURCE_DIR}/include/rpt/*.h")

add_library(lethe-rpt ${sources})

deal_ii_setup_target(lethe-rpt)
target_link_libraries(lethe-rpt lethe-core)
