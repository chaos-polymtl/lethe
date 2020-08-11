SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_nitsche_navier_stokes_23)
file(COPY noslip_23.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY taylor_couette_big_cylinder.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY taylor_couette_small_cylinder.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
DEAL_II_PICKUP_TESTS()

