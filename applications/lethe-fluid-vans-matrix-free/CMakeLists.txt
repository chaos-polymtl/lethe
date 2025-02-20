add_executable(lethe-fluid-vans-matrix-free fluid_dynamics_vans_matrix_free.cc)
deal_ii_setup_target(lethe-fluid-vans-matrix-free)
target_link_libraries(lethe-fluid-vans-matrix-free lethe-solvers lethe-fem-dem)
