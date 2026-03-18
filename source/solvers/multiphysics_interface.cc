// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "solvers/time_harmonic_maxwell.h"
#include <solvers/cahn_hilliard.h>
#include <solvers/cls.h>
#include <solvers/heat_transfer.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/tracer.h>

#include <deal.II/base/exceptions.h>
#define _unused(x) ((void)(x))

DeclException1(
  ThermalBuoyancyWithoutFluidDynamicsError,
  bool,
  << std::boolalpha << "Thermal buoyancy force is activated (" << arg1
  << "), while fluid dynamics is not activated (false)." << std::endl
  << "Thermal buoyancy force cannot be activated without activating fluid dynamics.");

DeclException1(
  ThermalBuoyancyWithoutHeatTransferError,
  bool,
  << std::boolalpha << "Thermal buoyancy force is activated (" << arg1
  << "), while heat transfer is not activated (false)." << std::endl
  << "Thermal buoyancy force cannot be activated without activating heat transfer.");

DeclException1(
  MicrowaveHeatingWithoutElectromagneticsError,
  bool,
  << std::boolalpha << "Microwave heating is activated (" << arg1
  << "), while electromagnetics is not activated (false)." << std::endl
  << "Microwave heating cannot be activated without activating electromagnetics.");

DeclException1(
  MarangoniWithoutFluidDynamicsError,
  bool,
  << std::boolalpha << "Marangoni effect is activated (" << arg1
  << "), while fluid dynamics is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating fluid dynamics.");

DeclException1(
  MarangoniWithoutHeatTransferError,
  bool,
  << std::boolalpha << "Marangoni effect is activated (" << arg1
  << "), while heat transfer is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating heat transfer.");

DeclException1(
  MarangoniWithoutSurfaceTensionForceError,
  bool,
  << std::boolalpha << "Marangoni effect is activated (" << arg1
  << "), while surface tension force is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating surface tension force.");

DeclException1(
  SurfaceTensionForceWithoutCLSError,
  bool,
  << std::boolalpha << "Surface tension force is activated (" << arg1
  << "), while CLS is not activated (false)." << std::endl
  << "Surface tension force cannot be activated without activating CLS.");

DeclException1(
  MarangoniWithoutCLSError,
  bool,
  << std::boolalpha << "Marangoni effect is activated (" << arg1
  << "), while CLS is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating CLS.");

DeclException1(
  InterfaceSharpeningWithoutCLSError,
  bool,
  << std::boolalpha << "Interface sharpening is activated (" << arg1
  << "), while CLS is not activated (false)." << std::endl
  << "Interface sharpening cannot be activated without activating CLS.");

DeclExceptionMsg(CahnHilliardWithHeatTransferError,
                 "Cahn-Hilliard and heat transfer are both activated. "
                 "This combination is not currently supported.");

DeclExceptionMsg(CahnHilliardWithThermalBuoyancyForceError,
                 "Cahn-Hilliard and thermal buoyancy force are both activated. "
                 "This combination is not currently supported.");

template <int dim>
MultiphysicsInterface<dim>::MultiphysicsInterface(
  const SimulationParameters<dim>                             &nsparam,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control,
  ConditionalOStream                &p_pcout)
  : multiphysics_parameters(nsparam.multiphysics)
  , simulation_control(p_simulation_control)
  , pcout(p_pcout)
{
  inspect_multiphysics_models_dependencies(nsparam);

  // Fluid dynamics is always considered active
  // since its DofHandler is required at all time by
  // the other physics. Consequently, disabling it only
  // prevents solving it, but not allocating it.
  {
    active_physics.push_back(PhysicsID::fluid_dynamics);
  }
  if (multiphysics_parameters.heat_transfer)
    {
      verbosity[PhysicsID::heat_transfer] =
        (nsparam.physics_solving_strategy.at(PhysicsID::heat_transfer)
             .verbosity != Parameters::Verbosity::quiet ||
         nsparam.linear_solver.at(PhysicsID::heat_transfer).verbosity !=
           Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;
      active_physics.push_back(PhysicsID::heat_transfer);
      physics[PhysicsID::heat_transfer] = std::make_shared<HeatTransfer<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }
  if (multiphysics_parameters.tracer)
    {
      verbosity[PhysicsID::tracer] =
        (nsparam.physics_solving_strategy.at(PhysicsID::tracer).verbosity !=
           Parameters::Verbosity::quiet ||
         nsparam.linear_solver.at(PhysicsID::tracer).verbosity !=
           Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;
      active_physics.push_back(PhysicsID::tracer);
      physics[PhysicsID::tracer] = std::make_shared<Tracer<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }
  if (multiphysics_parameters.CLS)
    {
      verbosity[PhysicsID::CLS] =
        (nsparam.physics_solving_strategy.at(PhysicsID::CLS).verbosity !=
           Parameters::Verbosity::quiet ||
         nsparam.linear_solver.at(PhysicsID::CLS).verbosity !=
           Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;
      active_physics.push_back(PhysicsID::CLS);
      physics[PhysicsID::CLS] = std::make_shared<ConservativeLevelSet<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }

  if (multiphysics_parameters.cahn_hilliard)
    {
      verbosity[PhysicsID::cahn_hilliard] =
        (nsparam.physics_solving_strategy.at(PhysicsID::cahn_hilliard)
             .verbosity != Parameters::Verbosity::quiet ||
         nsparam.linear_solver.at(PhysicsID::cahn_hilliard).verbosity !=
           Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;
      active_physics.push_back(PhysicsID::cahn_hilliard);
      physics[PhysicsID::cahn_hilliard] = std::make_shared<CahnHilliard<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }

  if (multiphysics_parameters.electromagnetics)
    {
      verbosity[PhysicsID::electromagnetics] =
        (nsparam.linear_solver.at(PhysicsID::electromagnetics).verbosity !=
         Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;

      active_physics.push_back(PhysicsID::electromagnetics);
      physics[PhysicsID::electromagnetics] =
        std::make_shared<TimeHarmonicMaxwell<dim>>(this,
                                                   nsparam,
                                                   p_triangulation,
                                                   p_simulation_control);
    }
}

template <int dim>
const GlobalVectorType &
MultiphysicsInterface<dim>::get_projected_phase_indicator_gradient_solution()
{
  // Throw error if CLS is not enabled
  AssertThrow((std::ranges::find(active_physics, PhysicsID::CLS) !=
               active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.cls_parameters.surface_tension_force.enable),
    ExcInternalError());

  return dynamic_cast<ConservativeLevelSet<dim> &>(*physics[PhysicsID::CLS])
    .get_projected_phase_indicator_gradient_solution();
}

template <int dim>
const GlobalVectorType &
MultiphysicsInterface<dim>::get_curvature_solution()
{
  // Throw error if CLS is not enabled
  AssertThrow((std::ranges::find(active_physics, PhysicsID::CLS) !=
               active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.cls_parameters.surface_tension_force.enable),
    ExcInternalError());

  return dynamic_cast<ConservativeLevelSet<dim> &>(*physics[PhysicsID::CLS])
    .get_curvature_solution();
}

template <int dim>
const DoFHandler<dim> &
MultiphysicsInterface<dim>::get_curvature_dof_handler()
{
  // Throw error if CLS is not enabled
  AssertThrow((std::ranges::find(active_physics, PhysicsID::CLS) !=
               active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.cls_parameters.surface_tension_force.enable),
    ExcInternalError());

  return dynamic_cast<ConservativeLevelSet<dim> &>(*physics[PhysicsID::CLS])
    .get_curvature_dof_handler();
}

template <int dim>
const DoFHandler<dim> &
MultiphysicsInterface<dim>::get_projected_phase_indicator_gradient_dof_handler()
{
  // Throw error if CLS is not enabled
  AssertThrow((std::ranges::find(active_physics, PhysicsID::CLS) !=
               active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.cls_parameters.surface_tension_force.enable),
    ExcInternalError());

  return dynamic_cast<ConservativeLevelSet<dim> &>(*physics[PhysicsID::CLS])
    .get_projected_phase_indicator_gradient_dof_handler();
}

template <int dim>
std::shared_ptr<Shape<dim>>
MultiphysicsInterface<dim>::get_immersed_solid_shape()
{
  // This shared pointer is also used to check if an immersed boundary solid
  // method is being used, depending on whether the pointer is assigned.
  return immersed_solid_shape;
}

template <int dim>
void
MultiphysicsInterface<dim>::set_immersed_solid_shape(
  const std::shared_ptr<Shape<dim>> &shape)
{
  immersed_solid_shape = shape;
}

template <int dim>
void
MultiphysicsInterface<dim>::inspect_multiphysics_models_dependencies(
  const SimulationParameters<dim> &nsparam)
{
  bool thermal_buoyancy_force_enabled =
    nsparam.multiphysics.thermal_buoyancy_force;
  bool heat_transfer_enabled     = nsparam.multiphysics.heat_transfer;
  bool microwave_heating_enabled = nsparam.multiphysics.microwave_heating;
  bool marangoni_effect_enabled =
    nsparam.multiphysics.cls_parameters.surface_tension_force
      .enable_marangoni_effect;
  bool surface_tension_force_enabled =
    nsparam.multiphysics.cls_parameters.surface_tension_force.enable;
  bool fluid_dynamics_enabled = nsparam.multiphysics.fluid_dynamics;
  bool interface_sharpening_enabled =
    nsparam.multiphysics.cls_parameters.reinitialization_method.sharpening
      .enable;
  bool CLS_enabled              = nsparam.multiphysics.CLS;
  bool cahn_hilliard_enabled    = nsparam.multiphysics.cahn_hilliard;
  bool electromagnetics_enabled = nsparam.multiphysics.electromagnetics;

  // To avoid getting unused parameter warning
  _unused(thermal_buoyancy_force_enabled && heat_transfer_enabled &&
          marangoni_effect_enabled && surface_tension_force_enabled &&
          fluid_dynamics_enabled && interface_sharpening_enabled &&
          CLS_enabled && cahn_hilliard_enabled && electromagnetics_enabled);

  // Dependence of thermal buoyancy force on fluid dynamics
  AssertThrow(!(thermal_buoyancy_force_enabled == true &&
                fluid_dynamics_enabled == false),
              ThermalBuoyancyWithoutFluidDynamicsError(
                thermal_buoyancy_force_enabled));

  // Dependence of thermal buoyancy force on heat transfer
  AssertThrow(
    !(thermal_buoyancy_force_enabled == true && heat_transfer_enabled == false),
    ThermalBuoyancyWithoutHeatTransferError(thermal_buoyancy_force_enabled));

  // Dependence of Marangoni effect on fluid dynamics
  AssertThrow(!(marangoni_effect_enabled == true &&
                fluid_dynamics_enabled == false),
              MarangoniWithoutFluidDynamicsError(marangoni_effect_enabled));

  // Dependence of Marangoni effect on heat transfer
  AssertThrow(!(marangoni_effect_enabled == true &&
                heat_transfer_enabled == false),
              MarangoniWithoutHeatTransferError(marangoni_effect_enabled));

  // Dependence of Marangoni effect on surface tension force
  AssertThrow(!(marangoni_effect_enabled == true &&
                surface_tension_force_enabled == false),
              MarangoniWithoutSurfaceTensionForceError(
                marangoni_effect_enabled));

  // Dependence of surface tension force on CLS
  AssertThrow(!(surface_tension_force_enabled == true && CLS_enabled == false),
              SurfaceTensionForceWithoutCLSError(
                surface_tension_force_enabled));

  // Dependence of Marangoni effect on CLS
  AssertThrow(!(marangoni_effect_enabled == true && CLS_enabled == false),
              MarangoniWithoutCLSError(marangoni_effect_enabled));

  // Dependence of interface sharpening on CLS
  AssertThrow(!(interface_sharpening_enabled == true && CLS_enabled == false),
              InterfaceSharpeningWithoutCLSError(interface_sharpening_enabled));

  // Cahn-Hilliard does not support heat transfer
  AssertThrow(!(cahn_hilliard_enabled == true && heat_transfer_enabled == true),
              CahnHilliardWithHeatTransferError());

  // Cahn-Hilliard does not support thermal buoyancy force
  AssertThrow(!(cahn_hilliard_enabled == true &&
                thermal_buoyancy_force_enabled == true),
              CahnHilliardWithThermalBuoyancyForceError());

  // Dependence of Microwave heating in heat transfer on electromagnetics
  AssertThrow(
    !(microwave_heating_enabled == true && electromagnetics_enabled == false),
    MicrowaveHeatingWithoutElectromagneticsError(microwave_heating_enabled));
}

template <int dim>
bool
MultiphysicsInterface<dim>::should_solve_electromagnetics() const
{
  const Parameters::TimeHarmonicMaxwell<dim> &thm_parameters =
    this->multiphysics_parameters.time_harmonic_maxwell_parameters;

  // For steady simulations, each outer iteration (e.g., mesh adaptation
  // cycles) should solve electromagnetics regardless of time-coupling settings.
  if (this->simulation_control->get_assembly_method() ==
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      return true;
    }

  // Always solve at the first step of the simulation (simulation start as 0 but
  // it is then incremented before solving the physics for the first time, so
  // the first time this function is called, the step number is 1).
  if (this->simulation_control->get_step_number() == 1)
    {
      return true;
    }
  else
    {
      switch (thm_parameters.time_coupling_method)
        {
          case Parameters::TimeCouplingMethod::none:
            return false;

            // Solve only if we are at a multiple of the specified iteration
            // frequency. We substract 1 since the step number starts at 1.
            if ((this->simulation_control->get_step_number() - 1) %
                  thm_parameters.coupling_iteration ==
                0)
              return true;
            else
              return false;
          case Parameters::TimeCouplingMethod::time:
            // Solve only if the current time has passed a multiple of
            // the specified time frequency since the last time the
            // electromagnetics were solved. This is done by comparing the floor
            // of the current time divided by the time coupling parameter to the
            // floor of the previous time divided by the time coupling
            // parameter. If they are different, it means we have passed a
            // multiple of the time coupling parameter and we should solve the
            // electromagnetics.
            if (std::floor(this->simulation_control->get_current_time() /
                           thm_parameters.coupling_time) >
                std::floor(this->simulation_control->get_previous_time() /
                           thm_parameters.coupling_time))
              return true;
            else
              return false;
          case Parameters::TimeCouplingMethod::threshold:
            AssertThrow(
              false,
              ExcMessage(
                "Time coupling method 'threshold' is not yet implemented "
                "for the time-harmonic Maxwell solver since the physical "
                "properties only support a constant model."));
            return true;
          default:
            AssertThrow(false, ExcMessage("Unknown time coupling method."));
            return false;
        }
    }
}

template class MultiphysicsInterface<2>;
template class MultiphysicsInterface<3>;
