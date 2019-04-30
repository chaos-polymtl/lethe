# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "cylinder")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(cylinder cylinder.cc)
DEAL_II_SETUP_TARGET(cylinder)
TARGET_LINK_LIBRARIES(cylinder lethe)

install(TARGETS cylinder RUNTIME DESTINATION bin)

#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
