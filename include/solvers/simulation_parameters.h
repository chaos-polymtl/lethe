/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef lethe_navier_stokes_solver_parameters_h
#define lethe_navier_stokes_solver_parameters_h

#include <core/boundary_conditions.h>
#include <core/dimensionality.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/solid_objects_parameters.h>

#include <solvers/analytical_solutions.h>
#include <solvers/initial_conditions.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/source_terms.h>

template <int dim>
class SimulationParameters
{
public:
  Parameters::Testing                               test;
  Parameters::LinearSolver                          linear_solver;
  Parameters::NonLinearSolver                       non_linear_solver;
  Parameters::MeshAdaptation                        mesh_adaptation;
  Parameters::Mesh                                  mesh;
  Parameters::Dimensionality                        dimensionality;
  std::shared_ptr<Parameters::MeshBoxRefinement>    mesh_box_refinement;
  std::shared_ptr<Parameters::Nitsche<dim>>         nitsche;
  Parameters::SimulationControl                     simulation_control;
  Parameters::Timer                                 timer;
  Parameters::FEM                                   fem_parameters;
  Parameters::Forces                                forces_parameters;
  std::shared_ptr<Parameters::Laser<dim>>           laser_parameters;
  Parameters::PostProcessing                        post_processing;
  Parameters::Restart                               restart_parameters;
  Parameters::Manifolds                             manifolds_parameters;
  BoundaryConditions::NSBoundaryConditions<dim>     boundary_conditions;
  BoundaryConditions::HTBoundaryConditions<dim>     boundary_conditions_ht;
  BoundaryConditions::TracerBoundaryConditions<dim> boundary_conditions_tracer;
  BoundaryConditions::VOFBoundaryConditions<dim>    boundary_conditions_vof;
  BoundaryConditions::CahnHilliardBoundaryConditions<dim>
                                      boundary_conditions_cahn_hilliard;
  Parameters::InitialConditions<dim> *initial_condition;
  AnalyticalSolutions::AnalyticalSolution<dim> *analytical_solution;
  SourceTerms::SourceTerm<dim> *                source_term;
  Parameters::VelocitySource                    velocity_sources;
  std::shared_ptr<Parameters::IBParticles<dim>> particlesParameters;
  Parameters::DynamicFlowControl                flow_control;
  Parameters::Multiphysics                      multiphysics;
  Parameters::Stabilization                     stabilization;

  PhysicalPropertiesManager physical_properties_manager;

  void
  declare(ParameterHandler &prm)
  {
    dimensionality.declare_parameters(prm);
    Parameters::SimulationControl::declare_parameters(prm);
    physical_properties.declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    nitsche = std::make_shared<Parameters::Nitsche<dim>>();
    nitsche->declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    boundary_conditions.declare_parameters(prm);
    boundary_conditions_ht.declare_parameters(prm);
    boundary_conditions_tracer.declare_parameters(prm);
    boundary_conditions_vof.declare_parameters(prm);
    boundary_conditions_cahn_hilliard.declare_parameters(prm);

    initial_condition = new Parameters::InitialConditions<dim>;
    initial_condition->declare_parameters(prm);

    Parameters::FEM::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Forces::declare_parameters(prm);
    laser_parameters = std::make_shared<Parameters::Laser<dim>>();
    laser_parameters->declare_parameters(prm);
    Parameters::MeshAdaptation::declare_parameters(prm);
    mesh_box_refinement = std::make_shared<Parameters::MeshBoxRefinement>();
    mesh_box_refinement->declare_parameters(prm);

    Parameters::NonLinearSolver::declare_parameters(prm);
    Parameters::LinearSolver::declare_parameters(prm);
    Parameters::PostProcessing::declare_parameters(prm);
    Parameters::DynamicFlowControl ::declare_parameters(prm);
    particlesParameters = std::make_shared<Parameters::IBParticles<dim>>();
    particlesParameters->declare_parameters(prm);
    manifolds_parameters.declare_parameters(prm);

    analytical_solution = new AnalyticalSolutions::AnalyticalSolution<dim>;
    analytical_solution->declare_parameters(prm);
    source_term = new SourceTerms::SourceTerm<dim>;
    source_term->declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);

    Parameters::VelocitySource::declare_parameters(prm);

    Parameters::Stabilization::declare_parameters(prm);

    multiphysics.declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    dimensionality.parse_parameters(prm);
    test.parse_parameters(prm);
    linear_solver.parse_parameters(prm);
    non_linear_solver.parse_parameters(prm);
    mesh_adaptation.parse_parameters(prm);
    mesh.parse_parameters(prm);
    mesh_box_refinement->parse_parameters(prm);
    nitsche->parse_parameters(prm);
    physical_properties.parse_parameters(prm, dimensionality);
    multiphysics.parse_parameters(prm);
    timer.parse_parameters(prm);
    fem_parameters.parse_parameters(prm);
    laser_parameters->parse_parameters(prm);
    forces_parameters.parse_parameters(prm);
    post_processing.parse_parameters(prm);
    flow_control.parse_parameters(prm);
    restart_parameters.parse_parameters(prm);
    boundary_conditions.parse_parameters(prm);
    boundary_conditions_ht.parse_parameters(prm);
    boundary_conditions_tracer.parse_parameters(prm);
    boundary_conditions_vof.parse_parameters(prm);
    boundary_conditions_cahn_hilliard.parse_parameters(prm);
    manifolds_parameters.parse_parameters(prm);
    initial_condition->parse_parameters(prm);
    analytical_solution->parse_parameters(prm);
    source_term->parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    velocity_sources.parse_parameters(prm);
    particlesParameters->parse_parameters(prm);
    multiphysics.parse_parameters(prm);
    stabilization.parse_parameters(prm);

    physical_properties_manager.initialize(physical_properties);


    // Check consistency of parameters parsed in different subsections
    if (multiphysics.VOF && physical_properties.number_of_fluids != 2)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n with VOF = true\n use: number of fluids = 2");
      }

    if (not(multiphysics.VOF) && post_processing.postprocessed_fluid ==
                                   Parameters::FluidIndicator::fluid1)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n when VOF = false"
          "\n use (default value): set postprocessed fluid = both"
          "\n or: set postprocessed fluid = fluid 0");
      }

    if (multiphysics.vof_parameters.sharpening.type ==
          Parameters::SharpeningType::adaptative &&
        not(multiphysics.vof_parameters.conservation.monitoring))
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n in subsection VOF, with sharpening type = adaptative\n "
          "use: monitoring = true");
      }
    if (multiphysics.vof_parameters.surface_tension_force.enable &&
        physical_properties.number_of_material_interactions == 0)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n in subsection VOF, with surface tension force enabled,\n but no material interactions specified in\n subsection physical properties\n"
          "use:\n"
          "  set number of material interactions = 1\n"
          "  subsection material interaction 0\n"
          "    set type = fluid-fluid\n"
          "    subsection fluid-fluid interaction\n"
          "      set first fluid id              = 0\n"
          "      set second fluid id             = 1\n"
          "      set surface tension model       = constant\n"
          "      set surface tension coefficient = $value_of_coefficient\n"
          "    end\n"
          "  end");
      }
    if (multiphysics.vof_parameters.surface_tension_force.enable)
      {
        bool error = true;
        for (unsigned int i = 0;
             i < physical_properties.number_of_material_interactions;
             ++i)
          {
            if (physical_properties.material_interactions[i]
                  .material_interaction_type ==
                Parameters::MaterialInteractions::MaterialInteractionsType::
                  fluid_fluid)
              error = false;
          }
        if (error)
          {
            throw std::logic_error(
              "Inconsistency in .prm!\n in subsection VOF, with surface tension force enabled,\n but no fluid-fluid material interactions specified in\n subsection physical properties\n"
              "use:\n"
              "  subsection material interaction $material_interaction_id\n"
              "    set type = fluid-fluid\n"
              "    subsection fluid-fluid interaction\n"
              "      set first fluid id              = 0\n"
              "      set second fluid id             = 1\n"
              "      set surface tension model       = constant\n"
              "      set surface tension coefficient = $value_of_coefficient\n"
              "    end\n"
              "  end");
          }
      }
  }

private:
  Parameters::PhysicalProperties physical_properties;
};

#endif
