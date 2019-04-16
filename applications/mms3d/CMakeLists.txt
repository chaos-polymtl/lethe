# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "mms2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )

ADD_EXECUTABLE(mms3d mms3d.cc)
DEAL_II_SETUP_TARGET(mms3d)
TARGET_LINK_LIBRARIES(mms3d lethe)
