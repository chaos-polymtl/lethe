set(TEST_TARGET dem_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY sliding_restart_files/sliding_restart.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")
file(COPY sliding_restart_files/sliding_restart.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/sliding_restart.${_build_type}")

file(COPY periodic_gmsh_files/pipe.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/insert_periodic_boundary_gmsh.${_build_type}")
file(COPY load_balancing_solid_object_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/load_balancing_solid_object.${_build_type}")
file(COPY moving_float_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/moving_float.${_build_type}")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(dem_3d/packing_in_box.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/packing_in_box.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/packing_in_ball.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/packing_in_ball.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/packing_in_cylinder.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/packing_in_ball_dynamic_contact.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/sliding_in_box.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(dem_3d/load_balancing_solid_object.mpirun=2.debug PROPERTIES TIMEOUT 3000)	
endif()
