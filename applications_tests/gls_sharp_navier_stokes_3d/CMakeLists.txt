SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_sharp_navier_stokes_3d)

DEAL_II_PICKUP_TESTS()

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY check_point.ib_particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY check_point.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY check_point.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY check_point.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY check_point.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY check_point.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY ib_force.00.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type})
file(COPY particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/load_particles_from_file_test.${_build_type})

	
if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
       	set_tests_properties(gls_sharp_navier_stokes_3d/steady_couette_sphere.debug PROPERTIES TIMEOUT 2000) 
       	set_tests_properties(gls_sharp_navier_stokes_3d/check_point.debug PROPERTIES TIMEOUT 2000) 
       	set_tests_properties(gls_sharp_navier_stokes_3d/pp_contact_test.debug PROPERTIES TIMEOUT 2000) 
       	set_tests_properties(gls_sharp_navier_stokes_3d/pp_lubrication_test.debug PROPERTIES TIMEOUT 2000) 
       	set_tests_properties(gls_sharp_navier_stokes_3d/coupled_moving_stokes.debug PROPERTIES TIMEOUT 2000) 
       	set_tests_properties(gls_sharp_navier_stokes_3d/pw_contact_test.debug PROPERTIES TIMEOUT 2000) 
ENDIF()

