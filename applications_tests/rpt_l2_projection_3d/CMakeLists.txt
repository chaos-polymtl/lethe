SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET rpt_l2_projection_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY rpt_l2_projection_files/positions.detector DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_l2_projection.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()


if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(rpt_l2_projection_3d/rpt_l2_projection.debug PROPERTIES TIMEOUT 1500)
ENDIF()