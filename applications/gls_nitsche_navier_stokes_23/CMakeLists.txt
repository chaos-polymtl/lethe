INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(gls_nitsche_navier_stokes_23 gls_nitsche_navier_stokes_23.cc)
DEAL_II_SETUP_TARGET(gls_nitsche_navier_stokes_23)
TARGET_LINK_LIBRARIES(gls_nitsche_navier_stokes_23 lethe-core lethe-solvers)

install(TARGETS gls_nitsche_navier_stokes_23 RUNTIME DESTINATION bin)

