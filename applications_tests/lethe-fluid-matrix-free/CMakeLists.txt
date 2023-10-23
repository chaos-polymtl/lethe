set(TEST_TARGET lethe-fluid-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if ( NOT ( "${DEAL_II_VERSION_MINOR}" STREQUAL 5))
    set_tests_properties(lethe-fluid-matrix-free/mms3d_fe1_gcmg.debug PROPERTIES TIMEOUT 800)
    set_tests_properties(lethe-fluid-matrix-free/mms3d_fe1_lsmg.debug PROPERTIES TIMEOUT 800)
  endif()
endif()