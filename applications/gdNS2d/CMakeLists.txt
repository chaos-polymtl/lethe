# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "gdNS2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(gdNS2d gdNS2d.cc)
DEAL_II_SETUP_TARGET(gdNS2d)
TARGET_LINK_LIBRARIES(gdNS2d lethe)

install(TARGETS gdNS2d RUNTIME DESTINATION bin)

