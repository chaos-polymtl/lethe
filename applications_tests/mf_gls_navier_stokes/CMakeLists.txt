<<<<<<< HEAD:applications_tests/mf_gls_navier_stokes/CMakeLists.txt
set(TEST_TARGET mf_gls_navier_stokes)
=======
set(TEST_TARGET mf_navier_stokes_2d)
>>>>>>> Rename all files:applications_tests/mf_navier_stokes_2d/CMakeLists.txt

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()
