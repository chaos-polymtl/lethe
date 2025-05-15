set(TEST_TARGET lethe-fluid-vans-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-vans-matrix-free/tanh_vans_A.debug PROPERTIES TIMEOUT 2000)
endif()