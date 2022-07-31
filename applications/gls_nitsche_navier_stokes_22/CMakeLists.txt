INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(gls_nitsche_navier_stokes_22 gls_nitsche_navier_stokes_22.cc)
DEAL_II_SETUP_TARGET(gls_nitsche_navier_stokes_22)
TARGET_LINK_LIBRARIES(gls_nitsche_navier_stokes_22 lethe-core lethe-solvers)

install(TARGETS gls_nitsche_navier_stokes_22 RUNTIME DESTINATION bin)

