set(TEST_TARGET gd_navier_stokes_2d)

file(COPY 2d_channel.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

deal_ii_pickup_tests()
