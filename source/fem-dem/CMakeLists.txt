add_library(lethe-fem-dem
  # Sources
  cfd_dem_coupling.cc
  cfd_dem_coupling_matrix_free.cc
  cfd_dem_simulation_parameters.cc
  fluid_dynamics_sharp.cc
  fluid_dynamics_vans.cc
  fluid_dynamics_vans_matrix_free.cc
  fluid_dynamics_vans_matrix_free_operators.cc
  ib_particles_dem.cc
  particle_projector.cc
  postprocessing_cfd_dem.cc
  vans_assemblers.cc
  # Headers
  ../../include/fem-dem/cfd_dem_coupling.h
  ../../include/fem-dem/cfd_dem_coupling_matrix_free.h
  ../../include/fem-dem/cfd_dem_simulation_parameters.h
  ../../include/fem-dem/fluid_dynamics_sharp.h
  ../../include/fem-dem/fluid_dynamics_vans.h
  ../../include/fem-dem/fluid_dynamics_vans_matrix_free.h
  ../../include/fem-dem/fluid_dynamics_vans_matrix_free_operators.h
  ../../include/fem-dem/ib_particles_dem.h
  ../../include/fem-dem/particle_projector.h
  ../../include/fem-dem/postprocessing_cfd_dem.h
  ../../include/fem-dem/vans_assemblers.h)

deal_ii_setup_target(lethe-fem-dem)
target_link_libraries(lethe-fem-dem lethe-core lethe-dem lethe-solvers)
