ADD_EXECUTABLE(cfd_dem_coupling_3d cfd_dem_coupling_3d.cc)
DEAL_II_SETUP_TARGET(cfd_dem_coupling_3d)
TARGET_LINK_LIBRARIES(cfd_dem_coupling_3d lethe-core lethe-solvers lethe-dem lethe-fem-dem)
