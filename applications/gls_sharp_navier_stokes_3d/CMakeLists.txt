INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(gls_sharp_navier_stokes_3d gls_sharp_navier_stokes_3d.cc)
DEAL_II_SETUP_TARGET(gls_sharp_navier_stokes_3d)
TARGET_LINK_LIBRARIES(gls_sharp_navier_stokes_3d lethe-core lethe-solvers lethe-dem lethe-fem-dem)

install(TARGETS gls_sharp_navier_stokes_3d RUNTIME DESTINATION bin)
