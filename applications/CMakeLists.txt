# top level CMakeLists.txt

ENABLE_TESTING()
ADD_SUBDIRECTORY(gls_navier_stokes_2d)
ADD_SUBDIRECTORY(gls_navier_stokes_3d)
ADD_SUBDIRECTORY(gd_navier_stokes_2d)
ADD_SUBDIRECTORY(gd_navier_stokes_3d)
ADD_SUBDIRECTORY(initial_conditions)
ADD_SUBDIRECTORY(navier_stokes_parameter_template)
ADD_SUBDIRECTORY(dem_3d)
ADD_SUBDIRECTORY(c_dem_3d)
ADD_SUBDIRECTORY(dem_parameter_template)


