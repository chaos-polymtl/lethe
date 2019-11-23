# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "NS_parameter_template")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(NS_parameter_template ns_parameter_template.cc)
DEAL_II_SETUP_TARGET(NS_parameter_template)
TARGET_LINK_LIBRARIES(NS_parameter_template lethe-core lethe-solvers)

install(TARGETS NS_parameter_template RUNTIME DESTINATION bin)

