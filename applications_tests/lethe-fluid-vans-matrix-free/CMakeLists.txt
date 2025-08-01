set(TEST_TARGET lethe-fluid-vans-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed_no_drag.${_build_type}/serial/")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-vans-matrix-free/tanh_vans_A.debug PROPERTIES TIMEOUT 2000)
endif()