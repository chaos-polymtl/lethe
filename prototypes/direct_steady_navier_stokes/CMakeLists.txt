##
#  CMake script for the step-6 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "direct_steady_navier_stokes")

ADD_EXECUTABLE(${TARGET} direct_steady_navier_stokes.cc
  boundaryconditions.h
  exactsolutions.h
  forcingfunctions.h)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET}) 
