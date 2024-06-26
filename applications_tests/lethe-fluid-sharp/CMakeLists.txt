IF(DEAL_II_WITH_OPENCASCADE)
add_subdirectory(opencascade)
ENDIF()
set(TEST_TARGET lethe-fluid-sharp)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY check_point_files/check_point.ib_particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.ib_particles.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/check_point.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")
file(COPY check_point_files/ib_force.00.dat DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point.${_build_type}/mpirun=1/")

file(COPY particles.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/load_particles_from_file_test.${_build_type}/mpirun=1/")
file(COPY rbf_test_shape.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/flow_around_rbf.${_build_type}/mpirun=1/")
file(COPY helix_composite_shape DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavatappi_composite_test.${_build_type}/mpirun=1/")
file(COPY spooky_composite_shape DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/spooky_composite_test.${_build_type}/mpirun=1/")
file(COPY rbf_sphere.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rbf_sphere_mms.${_build_type}/mpirun=1/")
file(COPY rbf_sphere_adaptive.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/rbf_sphere_mms_adaptive.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-sharp/steady_couette_sphere.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/check_point.mpirun=1.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/pp_contact_test.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/pp_lubrication_test.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/coupled_moving_stokes.mpirun=2.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/pw_contact_test.debug PROPERTIES TIMEOUT 2000)
  set_tests_properties(lethe-fluid-sharp/two_non_sphere_contact.debug PROPERTIES TIMEOUT 2000)
endif()
