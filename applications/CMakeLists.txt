foreach(_target IN ITEMS
    lethe-fluid-particles
    lethe-particles
    lethe-particles-ray-tracing
    lethe-fluid-block
    lethe-fluid
    lethe-fluid-nitsche
    lethe-fluid-sharp
    lethe-fluid-vans
    lethe-fluid-matrix-free
    lethe-fluid-vans-matrix-free
    lethe-fluid-particles-matrix-free
    lethe-parameter-check
    dem_parameter_template
    navier_stokes_parameter_template
    )
  add_subdirectory(${_target})
  install(TARGETS ${_target})
endforeach()
