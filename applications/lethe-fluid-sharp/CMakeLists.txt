add_executable(lethe-fluid-sharp fluid_dynamics_sharp.cc)
deal_ii_setup_target(lethe-fluid-sharp)
target_link_libraries(lethe-fluid-sharp lethe-fem-dem)
