IF(DEAL_II_WITH_OPENCASCADE)
set(TEST_TARGET lethe-fluid-sharp)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

file(COPY sphere_cad.step DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/static_stokes_step.${_build_type}/mpirun=1/")

deal_ii_pickup_tests()
ENDIF()
