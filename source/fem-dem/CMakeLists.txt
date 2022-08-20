add_library(lethe-fem-dem
  # Sources
  cfd_dem_coupling.cc
  cfd_dem_simulation_parameters.cc
  gls_sharp_navier_stokes.cc
  gls_vans.cc
  ib_particles_dem.cc
  parameters_cfd_dem.cc
  vans_assemblers.cc
  # Headers
  ../../include/fem-dem/cfd_dem_coupling.h
  ../../include/fem-dem/cfd_dem_simulation_parameters.h
  ../../include/fem-dem/gls_sharp_navier_stokes.h
  ../../include/fem-dem/gls_vans.h
  ../../include/fem-dem/ib_particles_dem.h
  ../../include/fem-dem/parameters_cfd_dem.h
  ../../include/fem-dem/vans_assemblers.h)

deal_ii_setup_target(lethe-fem-dem)
target_link_libraries(lethe-fem-dem lethe-core lethe-dem lethe-solvers)
