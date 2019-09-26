# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "glsNS2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
ADD_EXECUTABLE(glsNS2d glsNS2d.cc)
DEAL_II_SETUP_TARGET(glsNS2d)
TARGET_LINK_LIBRARIES(glsNS2d lethe-core lethe-solvers)

install(TARGETS glsNS2d RUNTIME DESTINATION bin)

