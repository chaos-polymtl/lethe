add_executable(dem dem.cc)
deal_ii_setup_target(dem)
target_link_libraries(dem lethe-core lethe-dem)
