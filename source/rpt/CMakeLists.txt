add_library(lethe-rpt
  # Sources
  compartmentalization.cc
  parameters_rpt.cc
  parameters_cp.cc
  particle_detector_interactions.cc
  particle_visualization.cc
  rpt.cc
  rpt_calculating_parameters.cc
  rpt_cell_reconstruction.cc
  rpt_fem_reconstruction.cc
  # Headers
  ../../include/rpt/compartmentalization.h
  ../../include/rpt/detector.h
  ../../include/rpt/parameters_rpt.h
  ../../include/rpt/parameters_cp.h
  ../../include/rpt/particle_detector_interactions.h
  ../../include/rpt/particle_visualization.h
  ../../include/rpt/radioactive_particle.h
  ../../include/rpt/rpt.h
  ../../include/rpt/rpt_calculating_parameters.h
  ../../include/rpt/rpt_cell_reconstruction.h
  ../../include/rpt/rpt_fem_reconstruction.h
  ../../include/rpt/rpt_utilities.h)

deal_ii_setup_target(lethe-rpt)
target_link_libraries(lethe-rpt lethe-core)
