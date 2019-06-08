# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "bdf")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(bdf_01 bdf_01.cc)
DEAL_II_SETUP_TARGET(bdf_01)
TARGET_LINK_LIBRARIES(bdf_01 lethe)
add_test( "BDF" bdf_01)

