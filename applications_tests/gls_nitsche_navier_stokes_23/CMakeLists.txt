set(TEST_TARGET gls_nitsche_navier_stokes_23)

file(COPY cylinder.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

deal_ii_pickup_tests()
