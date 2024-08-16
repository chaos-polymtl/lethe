file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/source/solvers/*.cc" "${CMAKE_SOURCE_DIR}/include/solvers/*.h")

add_library(lethe-solvers ${sources})

deal_ii_setup_target(lethe-solvers)
target_link_libraries(lethe-solvers lethe-core)
