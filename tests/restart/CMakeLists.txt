# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "restart_01")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(restart_01 restart_01.cc)
DEAL_II_SETUP_TARGET(restart_01)
TARGET_LINK_LIBRARIES(restart_01 lethe)
add_test( "Restart" restart_01)

