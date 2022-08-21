add_executable(gls_vans_3d gls_vans_3d.cc)
deal_ii_setup_target(gls_vans_3d)
target_link_libraries(gls_vans_3d lethe-fem-dem)
