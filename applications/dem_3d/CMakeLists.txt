add_executable(dem_3d dem_3d.cc)
deal_ii_setup_target(dem_3d)
target_link_libraries(dem_3d lethe-core lethe-dem)
