INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(dem_parameter_template dem_parameter_template.cc)
DEAL_II_SETUP_TARGET(dem_parameter_template)
TARGET_LINK_LIBRARIES(dem_parameter_template lethe-core lethe-solvers)

install(TARGETS dem_parameter_template RUNTIME DESTINATION bin)

