set(TEST_TARGET lethe-fluid-particles-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-particles-matrix-free/semi-implicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-particles-matrix-free/implicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-particles-matrix-free/explicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1200)
endif()
