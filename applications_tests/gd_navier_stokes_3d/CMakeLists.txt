SET(TEST_LIBRARIES lethe)
SET(TEST_TARGET gd_navier_stokes_3d)
DEAL_II_PICKUP_TESTS()


if (CMAKE_BUILD_TYPE STREQUAL "Debug" )
    if (${DEAL_II_VERSION_MAJOR} LESS 10)
       	set_tests_properties(gd_navier_stokes_3d/mms3d_gd.mpirun=2.debug PROPERTIES TIMEOUT 2000) 
    else()
       	set_tests_properties(gd_navier_stokes_3d/mms3d_gd.mpirun=2.threads=0.debug PROPERTIES TIMEOUT 2000) 
	ENDIF()
ENDIF()

