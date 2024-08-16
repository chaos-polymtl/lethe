file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/source/fem-dem/*.cc" "${CMAKE_SOURCE_DIR}/include/fem-dem/*.h")

add_library(lethe-fem-dem ${sources})

deal_ii_setup_target(lethe-fem-dem)
target_link_libraries(lethe-fem-dem lethe-core lethe-dem lethe-solvers)
