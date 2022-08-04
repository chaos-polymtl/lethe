SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gls_navier_stokes_3d)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY cylinder_unstructured.msh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY poiseuille_restart_files/poiseuille_restart.averagevelocities DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.flowcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.pvdhandler DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.simulationcontrol DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.triangulation DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.triangulation_fixed.data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart.triangulation.info DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0000.pvtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0005.pvtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0010.pvtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0015.pvtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0000.0000.vtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0005.0000.vtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0010.0000.vtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})
file(COPY poiseuille_restart_files/poiseuille_restart-output.0015.0000.vtu DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/poiseuille_restart.${_build_type})


DEAL_II_PICKUP_TESTS()

