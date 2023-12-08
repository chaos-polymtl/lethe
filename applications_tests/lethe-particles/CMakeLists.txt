set(TEST_TARGET lethe-particles)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY circle_restart_files/circle_restart.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")
file(COPY circle_restart_files/circle_restart.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${_build_type}")

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

# CMAKE instructions for disabled test due to boost compatibility issues
#file(COPY moving_receptacle_files/restart.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.solid_object.00.displacement DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.solid_object.00.dof DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")
#file(COPY moving_receptacle_files/restart.solid_object.00.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-particles/packing_in_box.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_box.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_cylinder.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball_dynamic_contact.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/sliding_in_box.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/load_balancing_solid_object.mpirun=2.debug PROPERTIES TIMEOUT 3000)	
endif()


