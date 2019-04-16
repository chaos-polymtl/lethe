# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "mms2d")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(mms2d mms2d.cc)
DEAL_II_SETUP_TARGET(mms2d)
TARGET_LINK_LIBRARIES(mms2d lethe)


#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
