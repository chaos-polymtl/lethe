SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_vans_2d)
DEAL_II_PICKUP_TESTS()

if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(gls_vans_2d/case3_vans.debug PROPERTIES TIMEOUT 1200) 
ENDIF()
