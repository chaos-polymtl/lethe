INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(cfd_dem_coupling_2d cfd_dem_coupling_2d.cc)
DEAL_II_SETUP_TARGET(cfd_dem_coupling_2d)
TARGET_LINK_LIBRARIES(cfd_dem_coupling_2d lethe-core lethe-solvers lethe-dem lethe-fem-dem)

install(TARGETS cfd_dem_coupling_2d RUNTIME DESTINATION bin)

