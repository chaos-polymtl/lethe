add_executable(lethe-fluid-nitsche nitsche_navier_stokes.cc)
deal_ii_setup_target(lethe-fluid-nitsche)
target_link_libraries(lethe-fluid-nitsche lethe-solvers)
