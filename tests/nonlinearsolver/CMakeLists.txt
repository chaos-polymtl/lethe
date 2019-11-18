# Set the name of the project and target:
SET(TARGET "nonlinearsolver")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include
  )
ADD_EXECUTABLE(nonlinearsolver nonlinearsolver.cc)
DEAL_II_SETUP_TARGET(nonlinearsolver)
TARGET_LINK_LIBRARIES(nonlinearsolver lethe-core)
add_test( "NONLINEARSOLVER" nonlinearsolver)