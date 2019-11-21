# Set the name of the project and target:
SET(TARGET "non_linear_solver_01")

INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include
  )
ADD_EXECUTABLE(non_linear_solver_01 nonlinearsolver_01.cc)
DEAL_II_SETUP_TARGET(non_linear_solver_01)
TARGET_LINK_LIBRARIES(non_linear_solver_01 lethe-core)
add_test( "non_linear_solver_01" non_linear_solver_01)
