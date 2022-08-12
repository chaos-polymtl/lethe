ADD_EXECUTABLE(gls_vans_2d gls_vans_2d.cc)
DEAL_II_SETUP_TARGET(gls_vans_2d)
TARGET_LINK_LIBRARIES(gls_vans_2d lethe-core lethe-solvers lethe-dem lethe-fem-dem)

install(TARGETS gls_vans_2d RUNTIME DESTINATION bin)

