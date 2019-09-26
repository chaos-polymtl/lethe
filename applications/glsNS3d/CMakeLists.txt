# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "glsNS3d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(glsNS3d glsNS3d.cc)
DEAL_II_SETUP_TARGET(glsNS3d)
TARGET_LINK_LIBRARIES(glsNS3d lethe-core lethe-solvers)

install(TARGETS glsNS3d RUNTIME DESTINATION bin)
