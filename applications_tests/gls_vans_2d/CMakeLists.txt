set(TEST_TARGET gls_vans_2d)
deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(gls_vans_2d/case3_vans.debug PROPERTIES TIMEOUT 1200)
endif()
