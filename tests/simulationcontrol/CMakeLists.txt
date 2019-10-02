# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "simulationcontrol_01")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include
  )
ADD_EXECUTABLE(simulationcontrol_01 simulationcontrol_01.cc)
DEAL_II_SETUP_TARGET(simulationcontrol_01)
TARGET_LINK_LIBRARIES(simulationcontrol_01 lethe-core)
add_test( "SimulationControl" simulationcontrol_01)

