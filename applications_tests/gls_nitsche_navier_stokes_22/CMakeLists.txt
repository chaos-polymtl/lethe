SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_nitsche_navier_stokes_22)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY rotating-square-restart/restart.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_fluid.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_fluid.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_fluid.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_fluid.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol.cells_vertices_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol.particles_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol.triangulation_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol.triangulation_00.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol.vertices_positions_00 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol_particles_00.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY rotating-square-restart/restart_sol_triangulation_00.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/rotating-square-restart.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()

