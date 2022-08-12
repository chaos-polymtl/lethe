INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
SET (TEST_LIBRARIES lethe-core lethe-solvers)
DEAL_II_PICKUP_TESTS()
