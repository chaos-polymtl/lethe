add_executable(gls_vans gls_vans.cc)
deal_ii_setup_target(gls_vans)
target_link_libraries(gls_vans lethe-fem-dem)
