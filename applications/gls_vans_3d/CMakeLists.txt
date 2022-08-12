ADD_EXECUTABLE(gls_vans_3d gls_vans_3d.cc)
DEAL_II_SETUP_TARGET(gls_vans_3d)
TARGET_LINK_LIBRARIES(gls_vans_3d lethe-core lethe-solvers lethe-dem lethe-fem-dem)

install(TARGETS gls_vans_3d RUNTIME DESTINATION bin)

