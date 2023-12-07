foreach(_target IN ITEMS
    lethe-fluid-particles
    lethe-particles
    lethe-fluid-block
    lethe-fluid
    lethe-fluid-nitsche
    lethe-fluid-sharp
    lethe-fluid-vans
    lethe-fluid-matrix-free
    lethe-parameter-check
    lethe-rpt-3d
    lethe-rpt-cell-reconstruction-3d
    lethe-rpt-fem-reconstruction-3d
    lethe-rpt-l2-projection-3d
    dem_parameter_template
    initial_conditions
    navier_stokes_parameter_template
    )
  add_subdirectory(${_target})
  install(TARGETS ${_target})
endforeach()
