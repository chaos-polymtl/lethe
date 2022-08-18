##
#  CMake script for the step-6 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "dg_heat_equation")

ADD_EXECUTABLE(${TARGET} dg_heat_equation.cc)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET})
