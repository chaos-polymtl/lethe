SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_vans_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY packed_bed_files/dem.particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})
file(COPY packed_bed_files/dem.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/packed_bed.${CMAKE_LOWER_BUILD_TYPE})


DEAL_II_PICKUP_TESTS()

