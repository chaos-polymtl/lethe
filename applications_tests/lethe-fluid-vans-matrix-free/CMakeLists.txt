set(TEST_TARGET lethe-fluid-vans-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()
