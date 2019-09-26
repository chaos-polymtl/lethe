# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "pvdhandler_01")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include
  )
ADD_EXECUTABLE(pvdhandler_01 pvdhandler_01.cc)
DEAL_II_SETUP_TARGET(pvdhandler_01)
TARGET_LINK_LIBRARIES(pvdhandler_01 lethe-core)
add_test( "PVDHandler" pvdhandler_01)

