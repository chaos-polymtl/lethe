SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET cfd_dem_coupling_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY particle_sedimentation_files/dem.particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})
file(COPY particle_sedimentation_files/dem.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/particle_sedimentation.${CMAKE_LOWER_BUILD_TYPE})


DEAL_II_PICKUP_TESTS()

