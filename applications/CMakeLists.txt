foreach(_target IN ITEMS
    cfd_dem_coupling_2d
    cfd_dem_coupling_3d
    compartmentalization_3d
    dem_2d
    dem_3d
    dem_parameter_template
    gd_navier_stokes_2d
    gd_navier_stokes_3d
    gls_navier_stokes_2d
    gls_navier_stokes_3d
    gls_nitsche_navier_stokes_22
    gls_nitsche_navier_stokes_23
    gls_nitsche_navier_stokes_33
    gls_sharp_navier_stokes_2d
    gls_sharp_navier_stokes_3d
    gls_vans_2d
    gls_vans_3d
    initial_conditions
    navier_stokes_parameter_template
    rpt_3d
    rpt_cell_reconstruction_3d
    rpt_fem_reconstruction_3d
    rpt_l2_projection_3d)
  add_subdirectory(${_target})
  install(TARGETS ${_target})
endforeach()
