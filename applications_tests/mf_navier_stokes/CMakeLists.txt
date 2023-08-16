set(TEST_TARGET mf_navier_stokes)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()
