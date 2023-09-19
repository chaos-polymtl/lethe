set(TEST_TARGET lethe-rpt-cell-reconstruction-3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY unknown_positions_files/counts.reconstruction DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/unknown_positions.${_build_type}")
file(COPY unknown_positions_files/positions.detector DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/unknown_positions.${_build_type}")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-rpt-cell-reconstruction-3d/unknown_positions.debug PROPERTIES TIMEOUT 1200)
endif()
