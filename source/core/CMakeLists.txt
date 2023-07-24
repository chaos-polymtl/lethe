add_library(lethe-core
  # Sources
  bdf.cc
  boundary_conditions.cc
  dem_properties.cc
  density_model.cc
  dimensionality.cc
  grids.cc
  ib_particle.cc
  ib_stencil.cc
  lethe_grid_tools.cc
  manifolds.cc
  mesh_controller.cc
  parameters.cc
  parameters_lagrangian.cc
  parameters_multiphysics.cc
  pvd_handler.cc
  rheological_model.cc
  sdirk.cc
  serial_solid.cc
  shape.cc
  shape_parsing.cc
  simulation_control.cc
  solid_base.cc
  solid_objects_parameters.cc
  solutions_output.cc
  specific_heat_model.cc
  surface_tension_model.cc
  thermal_conductivity_model.cc
  thermal_expansion_model.cc
  tracer_diffusivity_model.cc
  utilities.cc
  # Headers
  ../../include/core/auxiliary_math_functions.h
  ../../include/core/bdf.h
  ../../include/core/boundary_conditions.h
  ../../include/core/dem_properties.h
  ../../include/core/density_model.h
  ../../include/core/dimensionality.h
  ../../include/core/exceptions.h
  ../../include/core/grids.h
  ../../include/core/ib_particle.h
  ../../include/core/ib_stencil.h
  ../../include/core/inexact_newton_non_linear_solver.h
  ../../include/core/interface_property_model.h
  ../../include/core/kinsol_newton_non_linear_solver.h
  ../../include/core/lethe_grid_tools.h
  ../../include/core/manifolds.h
  ../../include/core/mesh_controller.h
  ../../include/core/multiphysics.h
  ../../include/core/newton_non_linear_solver.h
  ../../include/core/non_linear_solver.h
  ../../include/core/parameters.h
  ../../include/core/parameters_lagrangian.h
  ../../include/core/parameters_multiphysics.h
  ../../include/core/periodic_hills_grid.h
  ../../include/core/phase_change.h
  ../../include/core/physical_property_model.h
  ../../include/core/physics_solver.h
  ../../include/core/pvd_handler.h
  ../../include/core/rheological_model.h
  ../../include/core/sdirk.h
  ../../include/core/serial_solid.h
  ../../include/core/shape.h
  ../../include/core/shape_parsing.h
  ../../include/core/simulation_control.h
  ../../include/core/solid_base.h
  ../../include/core/solid_objects_parameters.h
  ../../include/core/solutions_output.h
  ../../include/core/specific_heat_model.h
  ../../include/core/surface_tension_model.h
  ../../include/core/tensors_and_points_dimension_manipulation.h
  ../../include/core/thermal_conductivity_model.h
  ../../include/core/thermal_expansion_model.h
  ../../include/core/time_integration_utilities.h
  ../../include/core/tracer_diffusivity_model.h
  ../../include/core/utilities.h)

deal_ii_setup_target(lethe-core)
target_include_directories(lethe-core PUBLIC ../../include)

# Don't add Boost in TARGET_LINK_LIBRARIES (here and in the other Lethe
# libraries) because it's added anyway by DEAL_II_SETUP_TARGET as part
# of deal.II's dependencies.
