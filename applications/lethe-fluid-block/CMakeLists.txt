add_executable(lethe-fluid-block gd_navier_stokes.cc)
deal_ii_setup_target(lethe-fluid-block)
target_link_libraries(lethe-fluid-block lethe-solvers)
