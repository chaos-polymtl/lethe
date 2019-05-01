# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "initial_conditions_L2")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(initial_conditions_L2 initial_conditions_L2.cc)
DEAL_II_SETUP_TARGET(initial_conditions_L2)
TARGET_LINK_LIBRARIES(initial_conditions_L2 lethe)

install(TARGETS initial_conditions_L2 RUNTIME DESTINATION bin)
