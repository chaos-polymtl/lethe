add_executable(lethe-particles dem.cc)
deal_ii_setup_target(lethe-particles)
target_link_libraries(lethe-particles lethe-core lethe-dem)
