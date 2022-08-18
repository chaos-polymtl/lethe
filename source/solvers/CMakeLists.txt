
# Set the name of the project and target:
SET(TARGET "lethe-solvers")

ADD_LIBRARY(${TARGET}
  # Sources
  analytical_solutions.cc
  auxiliary_physics.cc
  flow_control.cc
  gd_navier_stokes.cc
  gls_navier_stokes.cc
  gls_nitsche_navier_stokes.cc
  heat_transfer.cc
  heat_transfer_assemblers.cc
  heat_transfer_scratch_data.cc
  initial_conditions.cc
  multiphysics_interface.cc
  navier_stokes_assemblers.cc
  navier_stokes_base.cc
  navier_stokes_scratch_data.cc
  navier_stokes_vof_assemblers.cc
  physical_properties_manager.cc
  postprocessing_cfd.cc
  postprocessing_velocities.cc
  simulation_parameters.cc
  source_terms.cc
  tracer.cc
  tracer_assemblers.cc
  tracer_scratch_data.cc
  vof.cc
  vof_assemblers.cc
  vof_scratch_data.cc
  # Headers
  ../../include/solvers/analytical_solutions.h
  ../../include/solvers/auxiliary_physics.h
  ../../include/solvers/copy_data.h
  ../../include/solvers/flow_control.h
  ../../include/solvers/gd_navier_stokes.h
  ../../include/solvers/gls_navier_stokes.h
  ../../include/solvers/gls_nitsche_navier_stokes.h
  ../../include/solvers/heat_transfer.h
  ../../include/solvers/heat_transfer_assemblers.h
  ../../include/solvers/heat_transfer_scratch_data.h
  ../../include/solvers/initial_conditions.h
  ../../include/solvers/multiphysics_interface.h
  ../../include/solvers/navier_stokes_assemblers.h
  ../../include/solvers/navier_stokes_base.h
  ../../include/solvers/navier_stokes_scratch_data.h
  ../../include/solvers/navier_stokes_vof_assemblers.h
  ../../include/solvers/physical_properties_manager.h
  ../../include/solvers/post_processors.h
  ../../include/solvers/postprocessing_cfd.h
  ../../include/solvers/postprocessing_velocities.h
  ../../include/solvers/simulation_parameters.h
  ../../include/solvers/source_terms.h
  ../../include/solvers/stabilization.h
  ../../include/solvers/tracer.h
  ../../include/solvers/tracer_assemblers.h
  ../../include/solvers/tracer_scratch_data.h
  ../../include/solvers/vof.h
  ../../include/solvers/vof_assemblers.h
  ../../include/solvers/vof_scratch_data.h)

DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(lethe-solvers lethe-core)
install(TARGETS lethe-solvers ARCHIVE DESTINATION lib)

