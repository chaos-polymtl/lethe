SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET dem_2d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)

file(COPY circle_restart_files/circle_restart.particles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})
file(COPY circle_restart_files/circle_restart.triangulation_variable.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/circle_restart.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()

