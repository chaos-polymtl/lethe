set(TEST_TARGET lethe-fluid-particles-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
endif()
