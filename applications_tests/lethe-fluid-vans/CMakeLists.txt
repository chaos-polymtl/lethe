set(TEST_TARGET lethe-fluid-vans)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/qcm_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/difelice_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rong_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/dallavalle_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/kochhill_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/beetstra_packed_bed.${_build_type}/mpirun=1/")

file(COPY packed_bed_files/dem.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")
file(COPY packed_bed_files/dem_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gidaspow_packed_bed.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-vans/pcm_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/spm_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/qcm_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1600)
  set_tests_properties(lethe-fluid-vans/difelice_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/rong_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/dallavalle_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/kochhill_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/beetstra_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
  set_tests_properties(lethe-fluid-vans/gidaspow_packed_bed.mpirun=1.debug PROPERTIES TIMEOUT 1200)
endif()
