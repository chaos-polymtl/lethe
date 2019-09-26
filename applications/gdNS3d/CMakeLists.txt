# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "gdNS3d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )

ADD_EXECUTABLE(gdNS3d gdNS3d.cc)
DEAL_II_SETUP_TARGET(gdNS3d)
TARGET_LINK_LIBRARIES(gdNS3d lethe-core lethe-solvers)

install(TARGETS gdNS3d RUNTIME DESTINATION bin)
