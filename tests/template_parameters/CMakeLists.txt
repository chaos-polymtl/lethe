# Set the name of the project and target:
SET(TARGET "template_glsns")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include
  )
ADD_EXECUTABLE(template_glsns template_glsns.cc)
DEAL_II_SETUP_TARGET(template_glsns)
TARGET_LINK_LIBRARIES(template_glsns lethe-core lethe-solvers)
add_test( "template_GLSNS" template_glsns)

