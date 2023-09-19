add_executable(lethe-fluid-particles cfd_dem_coupling.cc)
deal_ii_setup_target(lethe-fluid-particles)
target_link_libraries(lethe-fluid-particles lethe-fem-dem)
