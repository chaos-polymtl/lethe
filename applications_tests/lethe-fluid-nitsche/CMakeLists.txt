set(TEST_TARGET lethe-fluid-nitsche)


string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY sphere.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

file(COPY two-bar-mixer-restart/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart.triangulation_variable.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_particles_00.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_particles_01.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_00.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_01.particles DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_00.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_00.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_00.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_01.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_01.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_01.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_triangulation_00.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")
file(COPY two-bar-mixer-restart/restart_solid_triangulation_01.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set_tests_properties(lethe-fluid-nitsche/noslip_22.debug PROPERTIES TIMEOUT 1200)
endif()

