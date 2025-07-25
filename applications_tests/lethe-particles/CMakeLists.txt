set(TEST_TARGET lethe-particles)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY restart_cell_weight_function_files/restart_cell_weight_function.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")
file(COPY restart_cell_weight_function_files/restart_cell_weight_function_1.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_cell_weight_function.${_build_type}/mpirun=2/")

file(COPY restart_circle_files/restart_circle.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")
file(COPY restart_circle_files/restart_circle_0.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_circle.${_build_type}/mpirun=1/")

file(COPY restart_sliding_files/restart_sliding.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")
file(COPY restart_sliding_files/restart_sliding_1.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_sliding.${_build_type}/mpirun=1/")

file(COPY periodic_gmsh_files/pipe.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/insert_periodic_boundary_gmsh.${_build_type}")
file(COPY load_balancing_solid_object_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/load_balancing_solid_object.${_build_type}")
file(COPY insert_file_3d_files/particles.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/insert_file_3d.${_build_type}")
file(COPY moving_solid_surface_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/moving_solid_surface_dmt.${_build_type}")
file(COPY moving_solid_surface_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/moving_solid_surface_hmlo.${_build_type}")
file(COPY moving_solid_surface_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/moving_solid_surface_jkr.${_build_type}")
file(COPY insert_and_remove_with_files_files/particles_00.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/insert_and_remove_with_files.${_build_type}")
file(COPY insert_and_remove_with_files_files/particles_01.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/insert_and_remove_with_files.${_build_type}")

file(COPY restart_moving_receptacle_files/restart.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.solid_object.00.displacement DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.solid_object.00.dof DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")
file(COPY restart_moving_receptacle_files/restart_1.solid_object.00.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_moving_receptacle.${_build_type}/mpirun=1/")

file(COPY restart_file_insertion_files/restart_file_insertion.checkpoint_controller DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/restart_file_insertion_1.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")
file(COPY restart_file_insertion_files/particles_01.input DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_file_insertion.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()

# Add the Utilis functions
include(../Utils.cmake)

# Get all .prm and .output files in the current directory
file(GLOB prm_files "${CMAKE_CURRENT_SOURCE_DIR}/*.prm")
file(GLOB output_files "${CMAKE_CURRENT_SOURCE_DIR}/*.output")

foreach (prm_file ${prm_files})

  get_depends_on("${prm_file}")
  if (_depends_on)

    message("Found a DEPENDS-ON: ${_depends_on}")
    # Name of the .prm file without the directory and the .prm extension.
    # e.g. lethe/lethe/application_test/lethe-particles/pp_jkr_equilibrium.prm --> pp_jkr_equilibrium
    get_filename_component(test_name ${prm_file} NAME_WLE)
    message("Test name without the path the the .prm extension: ${test_name}")

    # set matching_output_file  as an empty list
    set(matching_output_file  "")
    # Find all the files in the output_files list that have the same test_name
    foreach (output_file ${output_files})
      # If no match is found, pos will be equal to -1
      string(FIND "${output_file}" "${test_name}" pos)
      if (NOT pos EQUAL -1)
        # Append the output_file in the matching_output_file list
        list(APPEND matching_output_files "${output_file}")
      endif ()
    endforeach ()

  # Loop over the matching_output_file
    foreach (matching_output_file  ${matching_output_files})
      message("Output file name ${matching_output_file}")
      # Check if the mpirun=N is in the output string. This can be empty
      get_mpi_count("${matching_output_file}")

      message("Found mpirun string:${_mpirun_string}")

      set(full_test_name "${TEST_TARGET}/${test_name}${_mpirun_string}.${_build_type}")
      message("Full test name: ${full_test_name}")
      message("Second test name: ${_depends_on}")

      message("Setting dependency: ${full_test_name} DEPENDS ${_depends_on}")
      set_tests_properties("${full_test_name}" PROPERTIES DEPENDS "${_depends_on}")
    endforeach ()
  endif ()

endforeach ()




if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-particles/packing_in_box.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_box.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball.mpirun=1.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball.mpirun=2.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_cylinder.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/packing_in_ball_dynamic_contact.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/sliding_in_box.debug PROPERTIES TIMEOUT 3000)
  set_tests_properties(lethe-particles/load_balancing_solid_object.mpirun=2.debug PROPERTIES TIMEOUT 3000)	
  set_tests_properties(lethe-particles/insert_file_3d.mpirun=2.debug PROPERTIES TIMEOUT 3000)
endif()


