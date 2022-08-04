SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_vans_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY pcm_packed_bed_files/dem.particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})
file(COPY pcm_packed_bed_files/dem.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pcm_packed_bed.${_build_type})

file(COPY spm_packed_bed_files/dem.particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})
file(COPY spm_packed_bed_files/dem.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/spm_packed_bed.${_build_type})

DEAL_II_PICKUP_TESTS()

if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(gls_vans_3d/pcm_packed_bed.debug PROPERTIES TIMEOUT 1200) 
       	set_tests_properties(gls_vans_3d/spm_packed_bed.debug PROPERTIES TIMEOUT 1200) 
ENDIF()
