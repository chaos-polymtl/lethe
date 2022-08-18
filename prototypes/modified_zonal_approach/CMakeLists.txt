##
#  CMake script for the step-6 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "modified_zonal_approach")

ADD_EXECUTABLE(${TARGET} modified_zonal_approach.cc
  mockup_systems.h
  write_distances.h)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET} lethe-core lethe-solvers)
