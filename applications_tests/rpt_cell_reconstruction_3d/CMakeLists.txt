SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET rpt_cell_reconstruction_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY unknown_positions_files/counts.reconstruction DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/unknown_positions.${CMAKE_LOWER_BUILD_TYPE})
file(COPY unknown_positions_files/positions.detector DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/unknown_positions.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()

