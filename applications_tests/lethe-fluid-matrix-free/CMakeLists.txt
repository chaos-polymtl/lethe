set(TEST_TARGET lethe-fluid-matrix-free)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cylinder_kelly_estimator.${_build_type}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cylinder_kelly_estimator_lsmg.${_build_type}")
file(COPY cylinder_close_outlet.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cylinder_outlet_bc.${_build_type}")
file(COPY cylinder_structured.msh DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/cylinder_transient_avg.${_build_type}")


file(COPY tgv_restart_bdf1/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/enstrophy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/kinetic_energy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")
file(COPY tgv_restart_bdf1/L2Error_FD.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tgv_restart_bdf1.${_build_type}")


file(COPY turbulent_taylor_couette_restart/restart.pvdhandler DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/restart.simulationcontrol DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/restart.triangulation DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/restart.triangulation_fixed.data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/restart.triangulation.info DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/enstrophy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")
file(COPY turbulent_taylor_couette_restart/kinetic_energy.checkpoint DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/turbulent_taylor_couette_restart.${_build_type}")

deal_ii_pickup_tests()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
if ( NOT ( "${DEAL_II_VERSION_MINOR}" STREQUAL 5))
set_tests_properties(lethe-fluid-matrix-free/mms3d_fe1_gcmg.debug PROPERTIES TIMEOUT 800)
set_tests_properties(lethe-fluid-matrix-free/mms3d_fe1_lsmg.debug PROPERTIES TIMEOUT 800)
endif()
endif()
