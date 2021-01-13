SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_navier_stokes_2d)

string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_LOWER_BUILD_TYPE)


file(COPY 2d_channel.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY square.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY cylinder_structured.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY taylorCouette.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY mms-transient-conduction-restart_files/restart.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${CMAKE_LOWER_BUILD_TYPE})
file(COPY mms-transient-conduction-restart_files/restart.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${CMAKE_LOWER_BUILD_TYPE})
file(COPY mms-transient-conduction-restart_files/restart.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${CMAKE_LOWER_BUILD_TYPE})
file(COPY mms-transient-conduction-restart_files/restart.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${CMAKE_LOWER_BUILD_TYPE})
file(COPY mms-transient-conduction-restart_files/restart.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${CMAKE_LOWER_BUILD_TYPE})

DEAL_II_PICKUP_TESTS()

