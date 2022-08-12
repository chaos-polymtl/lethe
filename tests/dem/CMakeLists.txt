INCLUDE_DIRECTORIES(
  lethe
  ${CMAKE_SOURCE_DIR}/include/
  )
SET (TEST_LIBRARIES lethe-core lethe-dem)

DEAL_II_PICKUP_TESTS()


