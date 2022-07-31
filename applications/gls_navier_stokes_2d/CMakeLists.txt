INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(gls_navier_stokes_2d gls_navier_stokes_2d.cc)
DEAL_II_SETUP_TARGET(gls_navier_stokes_2d)
TARGET_LINK_LIBRARIES(gls_navier_stokes_2d lethe-core lethe-solvers)

install(TARGETS gls_navier_stokes_2d RUNTIME DESTINATION bin)

