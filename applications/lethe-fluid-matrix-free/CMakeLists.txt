add_executable(lethe-fluid-matrix-free fluid_dynamics_matrix_free.cc)
deal_ii_setup_target(lethe-fluid-matrix-free)
target_link_libraries(lethe-fluid-matrix-free lethe-solvers)
