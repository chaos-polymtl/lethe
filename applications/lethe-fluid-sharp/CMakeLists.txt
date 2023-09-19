add_executable(lethe-fluid-sharp gls_sharp_navier_stokes.cc)
deal_ii_setup_target(lethe-fluid-sharp)
target_link_libraries(lethe-fluid-sharp lethe-fem-dem)
