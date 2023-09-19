add_executable(lethe-fluid-matrix-free mf_navier_stokes.cc)
deal_ii_setup_target(lethe-fluid-matrix-free)
target_link_libraries(lethe-fluid-matrix-free lethe-solvers)
