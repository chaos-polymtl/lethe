SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_nitsche_navier_stokes_22)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY two-bar-mixer-restart/restart.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.particles_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol_particles_00.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.particles_01 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol_particles_01.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.triangulation_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.triangulation_00.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol_triangulation_00.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.triangulation_01 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol.triangulation_01.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY two-bar-mixer-restart/restart_sol_triangulation_01.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/two-bar-mixer-restart.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()

