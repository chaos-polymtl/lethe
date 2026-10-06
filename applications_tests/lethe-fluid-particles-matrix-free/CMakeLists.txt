set(TEST_TARGET lethe-fluid-particles-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY restart_particle_insertion_files/restart.insertion_object DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart_particles.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")
file(COPY restart_particle_insertion_files/restart.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/restart_particle_insertion.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-particles-matrix-free/semi-implicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1600)
  set_tests_properties(lethe-fluid-particles-matrix-free/implicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1600)
  set_tests_properties(lethe-fluid-particles-matrix-free/explicit-particle-sedimentation.debug PROPERTIES TIMEOUT 1600)
endif()
