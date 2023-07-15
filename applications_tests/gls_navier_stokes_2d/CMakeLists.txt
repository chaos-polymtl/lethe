set(TEST_TARGET gls_navier_stokes_2d)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY 2d_channel.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY square.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY taylorCouette.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY box.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

file(COPY cavity_ht_average_velocity_restart_files/restart.averagevelocities DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")
file(COPY cavity_ht_average_velocity_restart_files/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")
file(COPY cavity_ht_average_velocity_restart_files/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")
file(COPY cavity_ht_average_velocity_restart_files/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")
file(COPY cavity_ht_average_velocity_restart_files/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")
file(COPY cavity_ht_average_velocity_restart_files/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cavity_ht_average_velocity_restart.${_build_type}")

file(COPY mms-transient-conduction-restart_files/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${_build_type}")
file(COPY mms-transient-conduction-restart_files/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${_build_type}")
file(COPY mms-transient-conduction-restart_files/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${_build_type}")
file(COPY mms-transient-conduction-restart_files/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${_build_type}")
file(COPY mms-transient-conduction-restart_files/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/mms-transient-conduction-restart_gls.${_build_type}")

file(COPY taylor-green-vortex-restart_gls_bdf1/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/enstrophy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")
file(COPY taylor-green-vortex-restart_gls_bdf1/kinetic_energy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/taylor-green-vortex-restart_gls_bdf1.${_build_type}")

deal_ii_pickup_tests()
