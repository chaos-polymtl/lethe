add_executable(lethe-fluid-nitsche fluid_dynamics_nitsche.cc)
deal_ii_setup_target(lethe-fluid-nitsche)
target_link_libraries(lethe-fluid-nitsche lethe-solvers)
