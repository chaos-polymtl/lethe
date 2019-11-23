# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "gd_navier_stokes_2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(gd_navier_stokes_2d gd_navier_stokes_2d.cc)
DEAL_II_SETUP_TARGET(gd_navier_stokes_2d)
TARGET_LINK_LIBRARIES(gd_navier_stokes_2d lethe-core lethe-solvers)

install(TARGETS gd_navier_stokes_2d RUNTIME DESTINATION bin)

