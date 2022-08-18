##
#  CMake script for the step-6 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "directSteadyNS")

ADD_EXECUTABLE(${TARGET} directSteadyNS.cc
  boundaryconditions.h
  exactsolutions.h
  forcingfunctions.h)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET}) 
