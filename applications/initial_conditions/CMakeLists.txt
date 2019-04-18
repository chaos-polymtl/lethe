# Usually, you will not need to modify anything beyond this point...
#Method of manufactured solution 2D

# Set the name of the project and target:
SET(TARGET "initial_conditions")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/src/
  )
ADD_EXECUTABLE(initial_conditions initial_conditions.cc)
DEAL_II_SETUP_TARGET(initial_conditions)
TARGET_LINK_LIBRARIES(initial_conditions lethe)


#ENABLE_TESTING()
#ADD_SUBDIRECTORY(tests)
