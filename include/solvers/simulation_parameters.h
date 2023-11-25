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

#include <core/ale.h>
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
  std::map<PhysicsID, Parameters::LinearSolver>     linear_solver;
  std::map<PhysicsID, Parameters::NonLinearSolver>  non_linear_solver;
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
  SourceTerms::SourceTerm<dim>                  source_term;
  Parameters::VelocitySource                    velocity_sources;
  std::shared_ptr<Parameters::IBParticles<dim>> particlesParameters;
  Parameters::DynamicFlowControl                flow_control;
  Parameters::Multiphysics                      multiphysics;
  Parameters::Stabilization                     stabilization;
  Parameters::ALE<dim>                          ale;
  Parameters::Evaporation                       evaporation;



  PhysicalPropertiesManager physical_properties_manager;

  void
  declare(ParameterHandler             &prm,
          Parameters::SizeOfSubsections size_of_subsections)
  {
    prm.declare_entry("dimension",
                      "0",
                      Patterns::Integer(),
                      "Dimension of the problem");

    dimensionality.declare_parameters(prm);
    Parameters::SimulationControl::declare_parameters(prm);
    physical_properties.declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    nitsche = std::make_shared<Parameters::Nitsche<dim>>();
    nitsche->declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    boundary_conditions.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_ht.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_tracer.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_vof.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_cahn_hilliard.declare_parameters(
      prm, size_of_subsections.boundary_conditions);

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
    for (auto physics_name : physics_names)
      {
        Parameters::LinearSolver::declare_parameters(prm, physics_name);
        Parameters::NonLinearSolver::declare_parameters(prm, physics_name);
      }

    Parameters::PostProcessing::declare_parameters(prm);
    Parameters::DynamicFlowControl ::declare_parameters(prm);
    particlesParameters = std::make_shared<Parameters::IBParticles<dim>>();
    particlesParameters->declare_parameters(prm);
    manifolds_parameters.declare_parameters(prm);

    analytical_solution = new AnalyticalSolutions::AnalyticalSolution<dim>;
    analytical_solution->declare_parameters(prm);
    source_term.declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);

    Parameters::VelocitySource::declare_parameters(prm);

    Parameters::Stabilization::declare_parameters(prm);

    ale.declare_parameters(prm);

    Parameters::Evaporation::declare_parameters(prm);

    multiphysics.declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    dimensionality.parse_parameters(prm);
    test.parse_parameters(prm);

    for (auto physics_name : physics_names)
      {
        PhysicsID physics_id = get_physics_id(physics_name);
        linear_solver[physics_id].parse_parameters(prm, physics_name);
        non_linear_solver[physics_id].parse_parameters(prm, physics_name);
      }

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
    source_term.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    velocity_sources.parse_parameters(prm);
    particlesParameters->parse_parameters(prm);
    multiphysics.parse_parameters(prm);
    stabilization.parse_parameters(prm);
    ale.parse_parameters(prm);
    evaporation.parse_parameters(prm);

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
          Parameters::SharpeningType::adaptive &&
        not(multiphysics.vof_parameters.conservation.monitoring))
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n in subsection VOF, with sharpening type = adaptive\n "
          "use: monitoring = true");
      }

    // Interface physical property models consistency check
    if (multiphysics.vof_parameters.surface_tension_force.enable)
      {
        std::string constant_surface_tension_model(
          "    subsection fluid-fluid interaction\n"
          "      set first fluid id              = 0\n"
          "      set second fluid id             = 1\n"
          "      set surface tension model       = constant\n"
          "      set surface tension coefficient = $value_of_coefficient\n"
          "    end\n");

        std::string linear_surface_tension_model(
          "    subsection fluid-fluid interaction\n"
          "      set first fluid id                              = 0\n"
          "      set second fluid id                             = 1\n"
          "      set surface tension model                       = linear\n"
          "      set surface tension coefficient                 = $value_of_coefficient\n"
          "      set temperature-driven surface tension gradient = $value_of_gradient\n"
          "    end\n");

        std::string phase_change_surface_tension_model(
          "    subsection fluid-fluid interaction\n"
          "      set first fluid id                              = 0\n"
          "      set second fluid id                             = 1\n"
          "      set surface tension model                       = phase change\n"
          "      set surface tension coefficient                 = $value_of_coefficient\n"
          "      set temperature-driven surface tension gradient = $value_of_gradient\n"
          "      set solidus temperature                         = $value_of_solidus_temperature\n"
          "      set liquidus temperature                        = $value_of_liquidus_temperature\n"
          "    end\n");
        if (!multiphysics.vof_parameters.surface_tension_force
               .enable_marangoni_effect) // constant surface tension model
          {
            if (physical_properties.number_of_material_interactions == 0)
              {
                throw std::logic_error(
                  "Inconsistency in .prm!\n "
                  "In subsection VOF, with surface tension force enabled,\n "
                  "but no material interactions specified in\n "
                  "subsection physical properties.\n "
                  "Use:\n\n"
                  "  set number of material interactions = 1\n"
                  "  subsection material interaction 0\n"
                  "    set type = fluid-fluid\n" +
                  constant_surface_tension_model + "  end\n");
              }
            else if (multiphysics.VOF &&
                     multiphysics
                       .heat_transfer) // disabled Marangoni effect error
              {
                if (!is_constant_surface_tension_model(
                      physical_properties.material_interactions))
                  {
                    throw std::logic_error(
                      "Inconsistency in .prm!\n "
                      "In subsection multiphysics, VOF and heat transfer enabled,\n "
                      "and in subsection physical properties, a non-constant surface\n "
                      "tension model, but Marangoni effect disabled in subsection\n "
                      "surface tension force of subsection VOF. This is necessary to account\n "
                      "for Marangoni effect. In subsection VOF, use:\n\n "
                      "  subsection surface tension force\n"
                      "    set enable                  = true\n"
                      "    set enable marangoni effect = true\n"
                      "  end\n");
                  }
              }
            else
              {
                if (no_fluid_fluid_interaction_error(
                      physical_properties.material_interactions))
                  {
                    throw std::logic_error(
                      "Inconsistency in .prm!\n "
                      "in subsection VOF, surface tension force enabled,\n "
                      "but no fluid-fluid material interactions specified in \n "
                      "subsection physical properties\n "
                      "Use:\n\n"
                      "  subsection material interaction $material_interaction_id\n"
                      "    set type = fluid-fluid\n" +
                      constant_surface_tension_model + "  end\n");
                  }
              }
          }
        else // non-constant surface tension model
          {
            if (physical_properties.number_of_material_interactions == 0)
              {
                throw std::logic_error(
                  "Inconsistency in .prm!\n "
                  "In subsection VOF, marangoni effect enabled,\n "
                  "but no material interactions specified in subsection physical\n "
                  "properties. This is necessary to account for Marangoni \n "
                  "effect. In subsection physical properties, use:\n\n"
                  "  set number of material interactions = 1\n"
                  "  subsection material interaction 0\n"
                  "    set type = fluid-fluid\n" +
                  linear_surface_tension_model +
                  "\n"
                  "or: \n\n"
                  "  set number of material interactions = 1\n"
                  "  subsection material interaction 0\n"
                  "    set type = fluid-fluid\n" +
                  phase_change_surface_tension_model + "  end\n");
              }
            else
              {
                if (no_fluid_fluid_interaction_error(
                      physical_properties.material_interactions))
                  {
                    throw std::logic_error(
                      "Inconsistency in .prm!\n "
                      "In subsection VOF, Marangoni effect enabled,\n "
                      "but no fluid-fluid material interactions specified in subsection\n "
                      "physical properties. This is necessary to account for Marangoni\n "
                      "effect. In subsection physical properties, use:\n\n"
                      "  subsection material interaction $material_interaction_id\n"
                      "    set type = fluid-fluid\n" +
                      linear_surface_tension_model +
                      "\n"
                      "or:\n\n"
                      "  set number of material interactions = 1\n"
                      "  subsection material interaction 0\n"
                      "    set type = fluid-fluid\n" +
                      phase_change_surface_tension_model + "  end\n");
                  }
                if (is_constant_surface_tension_model(
                      physical_properties.material_interactions))
                  {
                    throw std::logic_error(
                      "Inconsistency in .prm!\n "
                      "In subsection VOF, Marangoni effect enabled,\n "
                      "but a constant surface tension model is specified in subsection\n "
                      "physical properties. This is necessary to account for Marangoni \n "
                      "effect. In subsection physical properties, use:\n\n"
                      "  subsection material interaction $material_interaction_id\n"
                      "    set type = fluid-fluid\n" +
                      linear_surface_tension_model +
                      "\n"
                      "or:\n\n"
                      "  set number of material interactions = 1\n"
                      "  subsection material interaction 0\n"
                      "    set type = fluid-fluid\n" +
                      phase_change_surface_tension_model + "  end\n");
                  }
              }
          }
      }

    if (multiphysics.cahn_hilliard && multiphysics.VOF)
      {
        throw std::runtime_error(
          "Cannot solve a multiphase problem using VOF and Cahn-Hilliard at the same time");
      }

    if (multiphysics.cahn_hilliard &&
        physical_properties.number_of_material_interactions == 0)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n "
          "cahn hilliard = true, \n "
          "but no fluid-fluid material interactions specified in\n "
          "subsection physical properties\n "
          "Use:\n\n"
          "  subsection material interaction $material_interaction_id\n"
          "    set type = fluid-fluid\n"
          "    subsection fluid-fluid interaction\n"
          "      set first fluid id                     = 0\n"
          "      set second fluid id                    = 1\n"
          "      set cahn hilliard mobility model       = constant\n"
          "      set cahn hilliard mobility coefficient = $value_of_coefficient\n"
          "    end\n"
          "  end\n");
      }
    if (laser_parameters->activate_laser &&
        laser_parameters->laser_type ==
          Parameters::Laser<dim>::LaserType::heat_flux_vof_interface &&
        !multiphysics.VOF)
      {
        throw std::logic_error(
          "At the moment, the laser surface heat flux is not implemented for 1 fluid simulations."
          "Please enable the VOF auxiliary physic in the 'multiphysics' subsection, \n"
          "specify a 2nd fluid in the 'physical properties' subsection,\n"
          "and define appropriate initial conditions in the 'initial conditions' subsection.");
      }
  }

  inline bool
  no_fluid_fluid_interaction_error(
    std::vector<Parameters::MaterialInteractions> &material_interactions)
  {
    for (const Parameters::MaterialInteractions &material_interaction :
         material_interactions)
      {
        if (material_interaction.material_interaction_type ==
            Parameters::MaterialInteractions::MaterialInteractionsType::
              fluid_fluid)
          {
            return false;
          }
      }
    return true;
  }

  inline bool
  is_constant_surface_tension_model(
    std::vector<Parameters::MaterialInteractions> &material_interactions)
  {
    for (const Parameters::MaterialInteractions &material_interaction :
         material_interactions)
      {
        if (material_interaction.surface_tension_model !=
            Parameters::MaterialInteractions::SurfaceTensionModel::constant)
          {
            return false;
          }
      }
    return true;
  }

private:
  Parameters::PhysicalProperties physical_properties;
  // names for the physics supported by Lethe
  std::vector<std::string> physics_names = {"fluid dynamics",
                                            "heat transfer",
                                            "tracer",
                                            "VOF",
                                            "cahn hilliard"};
};

#endif
