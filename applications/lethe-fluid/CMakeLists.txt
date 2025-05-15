add_executable(lethe-fluid fluid_dynamics_matrix_based.cc)
deal_ii_setup_target(lethe-fluid)
target_link_libraries(lethe-fluid lethe-solvers)
