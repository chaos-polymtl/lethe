add_executable(dem_2d dem_2d.cc)
deal_ii_setup_target(dem_2d)
target_link_libraries(dem_2d lethe-core lethe-dem)
