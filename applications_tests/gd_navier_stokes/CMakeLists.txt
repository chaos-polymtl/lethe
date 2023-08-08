set(TEST_TARGET gd_navier_stokes)

file(COPY 2d_channel.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(gd_navier_stokes_/mms3d_gd.mpirun=2.debug PROPERTIES TIMEOUT 2000)
endif()

deal_ii_pickup_tests()
