SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET lethe-rpt-fem-reconstruction-3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY rpt_fem_reconstruction_files/positions.detector DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE}/mpirun=1/)
file(COPY rpt_fem_reconstruction_files/flower_trajectory.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE}/mpirun=1/)

DEAL_II_PICKUP_TESTS()


if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(lethe-rpt-fem-reconstruction-3d/rpt_fem_reconstruction.mpirun=1.debug PROPERTIES TIMEOUT 2100)
ENDIF()