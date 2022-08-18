set(TEST_TARGET gls_nitsche_navier_stokes_33)

file(COPY sphere.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

deal_ii_pickup_tests()
