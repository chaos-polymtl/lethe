foreach(_target IN ITEMS
    cfd_dem_coupling
    dem
    dem_parameter_template
    gd_navier_stokes
    gls_navier_stokes
    gls_nitsche_navier_stokes_22
    gls_nitsche_navier_stokes_23
    gls_nitsche_navier_stokes_33
    gls_sharp_navier_stokes
    gls_vans
    mf_gls_navier_stokes
    initial_conditions
    navier_stokes_parameter_template
    rpt_3d
    rpt_cell_reconstruction_3d
    rpt_fem_reconstruction_3d
    rpt_l2_projection_3d)
  add_subdirectory(${_target})
  install(TARGETS ${_target})
endforeach()
