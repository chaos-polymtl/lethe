set(TEST_TARGET lethe-fluid-sharp-filter)

string(TOLOWER ${CMAKE_BUILD_TYPE} _build_type)

# The restart test filters the checkpoint written by the check_point test of
# lethe-fluid-sharp.
foreach(_file IN ITEMS
    check_point.ib_particles
    check_point.ib_particles.pvdhandler
    check_point.pvdhandler
    check_point.simulationcontrol
    check_point.triangulation
    check_point.triangulation.info
    check_point.triangulation_fixed.data
    ib_force.00.dat)
  file(COPY ../lethe-fluid-sharp/check_point_files/${_file} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/check_point_restart.${_build_type}/mpirun=1/")
endforeach()

deal_ii_pickup_tests()
