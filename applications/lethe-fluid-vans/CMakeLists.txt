add_executable(lethe-fluid-vans fluid_dynamics_vans.cc)
deal_ii_setup_target(lethe-fluid-vans)
target_link_libraries(lethe-fluid-vans lethe-fem-dem)
