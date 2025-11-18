add_executable(lethe-fluid-particles-matrix-free lethe_fluid_particles_matrix_free.cc)
deal_ii_setup_target(lethe-fluid-particles-matrix-free)
target_link_libraries(lethe-fluid-particles-matrix-free lethe-fem-dem)
