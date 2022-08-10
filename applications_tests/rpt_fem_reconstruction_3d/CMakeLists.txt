SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET rpt_fem_reconstruction_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY rpt_fem_reconstruction_files/positions.detector DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/flower_trajectory.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_dof_handler.dof DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_dof_handler.dof DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector00.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector01.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector02.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector03.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector04.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rpt_fem_reconstruction_files/temp_nodal_counts_detector05.counts DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rpt_fem_reconstruction.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()


if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(rpt_fem_reconstruction_3d/rpt_fem_reconstruction.debug PROPERTIES TIMEOUT 300)
ENDIF()