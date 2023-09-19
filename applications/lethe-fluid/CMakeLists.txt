add_executable(lethe-fluid gls_navier_stokes.cc)
deal_ii_setup_target(lethe-fluid)
target_link_libraries(lethe-fluid lethe-solvers)
