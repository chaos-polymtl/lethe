# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "navier_stokes_parameter_template")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(navier_stokes_parameter_template navier_stokes_parameter_template.cc)
DEAL_II_SETUP_TARGET(navier_stokes_parameter_template)
TARGET_LINK_LIBRARIES(navier_stokes_parameter_template lethe-core lethe-solvers)

install(TARGETS navier_stokes_parameter_template RUNTIME DESTINATION bin)

