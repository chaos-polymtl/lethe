// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_simulation_parameters_h
#define lethe_simulation_parameters_h

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
#include <solvers/tracer_drift_velocity.h>
#include <solvers/cls_subequations.h>

template <int dim>
class SimulationParameters
{
public:
  Parameters::Testing                              test;
  std::map<PhysicsID, Parameters::LinearSolver>    linear_solver;
  std::map<PhysicsID, Parameters::NonLinearSolver> physics_solving_strategy;
  std::map<CLSSubequationsID, Parameters::LinearSolver>
    cls_subequations_linear_solvers;
  std::map<CLSSubequationsID, Parameters::NonLinearSolver>
                             cls_subequations_non_linear_solvers;
  Parameters::MeshAdaptation mesh_adaptation;
  Parameters::Mesh           mesh;
  Parameters::Dimensionality dimensionality;
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
  BoundaryConditions::CLSBoundaryConditions<dim>    boundary_conditions_cls;
  BoundaryConditions::CahnHilliardBoundaryConditions<dim>
    boundary_conditions_cahn_hilliard;
  BoundaryConditions::TimeHarmonicMaxwellBoundaryConditions<dim>
    boundary_conditions_time_harmonic_electromagnetics;
  Parameters::InitialConditions<dim>           *initial_condition;
  AnalyticalSolutions::AnalyticalSolution<dim> *analytical_solution;
  SourceTerms::SourceTerm<dim>                  source_term;
  Parameters::VelocitySource                    velocity_sources;
  std::shared_ptr<Parameters::IBParticles<dim>> particlesParameters;
  Parameters::DynamicFlowControl                flow_control;
  Parameters::Multiphysics<dim>                 multiphysics;
  Parameters::ConstrainSolidDomain<dim>         constrain_solid_domain;
  Parameters::Stabilization                     stabilization;
  Parameters::ALE<dim>                          ale;
  Parameters::Evaporation                       evaporation;
  Parameters::TracerDriftVelocity<dim>          tracer_drift_velocity;
  Parameters::Mortar<dim>                       mortar_parameters;


  PhysicalPropertiesManager physical_properties_manager;

  void
  declare(ParameterHandler             &prm,
          Parameters::SizeOfSubsections size_of_subsections)
  {
    prm.declare_entry("dimension",
                      "0",
                      Patterns::Integer(),
                      "Dimension of the problem");

    prm.declare_entry("print parameters",
                      "none",
                      Patterns::Selection("none|only changed|all"),
                      "Print all the parameters, or only"
                      "the changed parameters or none");

    prm.declare_entry(
      "comment message",
      "",
      Patterns::Anything(),
      "Print a comment at the beginning of the console output.");

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
    boundary_conditions_cls.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_cahn_hilliard.declare_parameters(
      prm, size_of_subsections.boundary_conditions);
    boundary_conditions_time_harmonic_electromagnetics.declare_parameters(
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
    for (auto physics_name : nonlinear_physics_names)
      {
        Parameters::LinearSolver::declare_parameters(prm, physics_name);
        Parameters::NonLinearSolver::declare_parameters(prm, physics_name);
      }
    for (auto physics_name : linear_physics_names)
      {
        Parameters::LinearSolver::declare_parameters(prm, physics_name);
      }
    for (const auto &cls_subequation_name : cls_subequations_names)
      {
        Parameters::LinearSolver::declare_parameters(prm, cls_subequation_name);
        Parameters::NonLinearSolver::declare_parameters(prm,
                                                        cls_subequation_name);
      }

    Parameters::PostProcessing::declare_parameters(prm);
    Parameters::DynamicFlowControl ::declare_parameters(prm);
    particlesParameters = std::make_shared<Parameters::IBParticles<dim>>();
    particlesParameters->declare_parameters(prm);
    manifolds_parameters.declare_parameters(prm, size_of_subsections.manifolds);

    analytical_solution = new AnalyticalSolutions::AnalyticalSolution<dim>;
    analytical_solution->declare_parameters(prm);
    source_term.declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);

    Parameters::VelocitySource::declare_parameters(prm);

    constrain_solid_domain.declare_parameters(prm, 2);

    Parameters::Stabilization::declare_parameters(prm);

    ale.declare_parameters(prm);

    Parameters::Evaporation::declare_parameters(prm);

    multiphysics.declare_parameters(prm);

    tracer_drift_velocity.declare_parameters(prm);

    mortar_parameters.declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    dimensionality.parse_parameters(prm);
    test.parse_parameters(prm);

    for (auto physics_name : nonlinear_physics_names)
      {
        PhysicsID physics_id = get_physics_id(physics_name);
        linear_solver[physics_id].parse_parameters(prm, physics_name);
        physics_solving_strategy[physics_id].parse_parameters(prm,
                                                              physics_name);
      }
    for (auto physics_name : linear_physics_names)
      {
        PhysicsID physics_id = get_physics_id(physics_name);
        linear_solver[physics_id].parse_parameters(prm, physics_name);
      }
    for (const auto &cls_subequation_name : cls_subequations_names)
      {
        CLSSubequationsID cls_subequations_id =
          get_cls_subequation_id(cls_subequation_name);
        cls_subequations_linear_solvers[cls_subequations_id].parse_parameters(
          prm, cls_subequation_name);
        cls_subequations_non_linear_solvers[cls_subequations_id]
          .parse_parameters(prm, cls_subequation_name);
      }

    mesh_adaptation.parse_parameters(prm);
    mesh.parse_parameters(prm);
    mesh_box_refinement->parse_parameters(prm);
    nitsche->parse_parameters(prm);
    physical_properties.parse_parameters(prm, dimensionality);
    multiphysics.parse_parameters(prm, dimensionality);
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
    boundary_conditions_cls.parse_parameters(prm);
    boundary_conditions_cahn_hilliard.parse_parameters(prm);
    boundary_conditions_time_harmonic_electromagnetics.parse_parameters(prm);
    manifolds_parameters.parse_parameters(prm);
    initial_condition->parse_parameters(prm);
    analytical_solution->parse_parameters(prm);
    source_term.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    velocity_sources.parse_parameters(prm);
    particlesParameters->parse_parameters(prm);
    constrain_solid_domain.parse_parameters(prm);
    stabilization.parse_parameters(prm);
    ale.parse_parameters(prm);
    evaporation.parse_parameters(prm);
    tracer_drift_velocity.parse_parameters(prm);
    mortar_parameters.parse_parameters(prm);

    physical_properties_manager.initialize(physical_properties);


    // Check consistency of parameters parsed in different subsections
    if (multiphysics.CLS && physical_properties.number_of_fluids != 2)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n with CLS = true\n use: number of fluids = 2");
      }

    if (not(multiphysics.CLS) && post_processing.postprocessed_fluid ==
                                   Parameters::FluidIndicator::fluid1)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n when CLS = false"
          "\n use (default value): set postprocessed fluid = both"
          "\n or: set postprocessed fluid = fluid 0");
      }

    if (physical_properties.number_of_fluids == 2 &&
        (!multiphysics.CLS && !multiphysics.cahn_hilliard))
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n "
          "Number of fluids in 'physical properties' was set to 2,\n "
          "but neither 'cls' or 'cahn hilliard' is enabled in the 'multiphysics'.\n ");
      }

    // Interface physical property models consistency check
    if (multiphysics.cls_parameters.surface_tension_force.enable)
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
        if (!multiphysics.cls_parameters.surface_tension_force
               .enable_marangoni_effect) // constant surface tension model
          {
            if (physical_properties.number_of_material_interactions == 0)
              {
                throw std::logic_error(
                  "Inconsistency in .prm!\n "
                  "In subsection CLS, with surface tension force enabled,\n "
                  "but no material interactions specified in\n "
                  "subsection physical properties.\n "
                  "Use:\n\n"
                  "  set number of material interactions = 1\n"
                  "  subsection material interaction 0\n"
                  "    set type = fluid-fluid\n" +
                  constant_surface_tension_model + "  end\n");
              }
            else if (multiphysics.CLS &&
                     multiphysics
                       .heat_transfer) // disabled Marangoni effect error
              {
                if (!is_constant_surface_tension_model(
                      physical_properties.material_interactions))
                  {
                    throw std::logic_error(
                      "Inconsistency in .prm!\n "
                      "In subsection multiphysics, 'cls' and 'heat transfer' enabled,\n "
                      "and in subsection physical properties, a non-constant surface\n "
                      "tension model, but Marangoni effect disabled in subsection\n "
                      "surface tension force of subsection CLS. This is necessary to account\n "
                      "for Marangoni effect. In subsection CLS, use:\n\n "
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
                      "in subsection CLS, surface tension force enabled,\n "
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
                  "In subsection CLS, marangoni effect enabled,\n "
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
                      "In subsection CLS, Marangoni effect enabled,\n "
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
                      "In subsection CLS, Marangoni effect enabled,\n "
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

    if (multiphysics.cahn_hilliard && multiphysics.CLS)
      {
        throw std::runtime_error(
          "Cannot solve a multiphase problem using CLS and Cahn-Hilliard at the same time");
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
        (laser_parameters->laser_type ==
           Parameters::Laser<
             dim>::LaserType::gaussian_heat_flux_cls_interface ||
         laser_parameters->laser_type ==
           Parameters::Laser<
             dim>::LaserType::uniform_heat_flux_cls_interface) &&
        !multiphysics.CLS)
      {
        throw std::logic_error(
          "At the moment, the laser surface heat flux is not implemented for 1 fluid simulations."
          "Please enable the CLS auxiliary physic in the 'multiphysics' subsection, \n"
          "specify a 2nd fluid in the 'physical properties' subsection,\n"
          "and define appropriate initial conditions in the 'initial conditions' subsection.");
      }

    if (!multiphysics.heat_transfer && constrain_solid_domain.enable)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n "
          "The apply constraints on a solid domain feature is enabled, however\n "
          "'heat transfer' was not set to 'true' in the 'multiphysics' subsection.\n ");
      }

    if (constrain_solid_domain.enable &&
        physical_properties_manager.get_number_of_fluids() <
          constrain_solid_domain.number_of_constraints)
      {
        std::string n_constraints =
          Utilities::to_string(constrain_solid_domain.number_of_constraints);
        std::string n_fluids = Utilities::to_string(
          physical_properties_manager.get_number_of_fluids());
        throw std::logic_error(
          "Inconsistency in .prm!\n "
          "The number of constraints (" +
          n_constraints + ") is greater than the number of fluids (" +
          n_fluids +
          ").\n "
          "Only 1 constraint per fluid can be declared.\n ");
      }

    if (constrain_solid_domain.enable && multiphysics.cahn_hilliard)
      {
        throw std::logic_error(
          "Inconsistency in .prm!\n "
          "The current implementation for constraining solid domains with\n "
          "temperature is not implemented for Cahn Hilliard simulations.\n ");
      }

    if (multiphysics.thermal_buoyancy_force)
      {
        for (const auto &fluid : physical_properties.fluids)
          {
            AssertThrow(
              fluid.density_model !=
                Parameters::Material::DensityModel::isothermal_ideal_gas,
              ExcMessage(
                "Inconsistency in .prm!\n "
                "Thermal buoyancy force is enabled while the isothermal ideal gas\n "
                "density model is used. The isothermal ideal gas density model assumes\n "
                "thermal independence of the density, which is incompatible with the\n "
                "thermal buoyancy force."));
          }
      }

    if (simulation_control.adapt_with_capillary_time_step_ratio)
      AssertThrow(
        (multiphysics.cls_parameters.surface_tension_force.enable &&
         multiphysics.CLS),
        ExcMessage(
          "The current implementation only allows the capillary time-step constraint \n "
          "to be respected for CLS multiphase flows with surface tension.\n "));
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
  std::vector<std::string> nonlinear_physics_names = {"fluid dynamics",
                                                      "heat transfer",
                                                      "tracer",
                                                      "CLS",
                                                      "cahn hilliard",
                                                      "void fraction"};
  std::vector<std::string> linear_physics_names    = {"electromagnetics"};
  // Names of subequations within CLS that inherits from PhysicsSolver
  std::vector<std::string> cls_subequations_names = {
    "CLS PDE-based interface reinitialization"};
};

#endif
