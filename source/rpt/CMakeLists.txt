
# Set the name of the project and target:
SET(TARGET "lethe-rpt")

ADD_LIBRARY(${TARGET}
  # Sources
  parameters_rpt.cc
  particle_detector_interactions.cc
  particle_visualization.cc
  rpt.cc
  rpt_calculating_parameters.cc
  rpt_cell_reconstruction.cc
  # Headers
  ../../include/rpt/detector.h
  ../../include/rpt/parameters_rpt.h
  ../../include/rpt/particle_detector_interactions.h
  ../../include/rpt/particle_visualization.h
  ../../include/rpt/radioactive_particle.h
  ../../include/rpt/rpt.h
  ../../include/rpt/rpt_calculating_parameters.h
  ../../include/rpt/rpt_cell_reconstruction.h
  ../../include/rpt/rpt_utilities.h)

DEAL_II_SETUP_TARGET(${TARGET})
TARGET_LINK_LIBRARIES(lethe-rpt lethe-core)
install(TARGETS lethe-rpt ARCHIVE DESTINATION lib)

