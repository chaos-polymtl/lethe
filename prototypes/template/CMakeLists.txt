##
#  CMake script for the step-6 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "protype_template")

ADD_EXECUTABLE(${TARGET} prototype_template.cc)
DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(${TARGET} lethe-core lethe-solvers)
