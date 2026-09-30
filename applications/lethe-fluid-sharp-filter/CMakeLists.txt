add_executable(lethe-fluid-sharp-filter lethe_fluid_sharp_filter.cc)
deal_ii_setup_target(lethe-fluid-sharp-filter)
target_link_libraries(lethe-fluid-sharp-filter lethe-fem-dem)
