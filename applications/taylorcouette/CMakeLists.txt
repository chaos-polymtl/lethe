# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "taylorcouette")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(taylorcouette taylorcouette.cc)
DEAL_II_SETUP_TARGET(taylorcouette)
TARGET_LINK_LIBRARIES(taylorcouette lethe)


#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
