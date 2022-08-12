add_executable(gls_vans_2d gls_vans_2d.cc)
deal_ii_setup_target(gls_vans_2d)
target_link_libraries(gls_vans_2d lethe-fem-dem)
