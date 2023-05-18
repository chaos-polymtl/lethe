#include <solvers/heat_transfer.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/tracer.h>
#include <solvers/vof.h>

#include <deal.II/base/exceptions.h>
#define _unused(x) ((void)(x))

DeclException1(
  BuoyancyWithoutFluidDynamicsError,
  bool,
  << "Buoyancy force is activated (" << arg1
  << "), while fluid dynamics is not activated (false)." << std::endl
  << "Buoyancy force cannot be activated without activating fluid dynamics.");

DeclException1(
  BuoyancyWithoutHeatTransferError,
  bool,
  << "Buoyancy force is activated (" << arg1
  << "), while heat transfer is not activated (false)." << std::endl
  << "Buoyancy force cannot be activated without activating heat transfer.");

DeclException1(
  MarangoniWithoutFluidDynamicsError,
  bool,
  << "Marangoni effect is activated (" << arg1
  << "), while fluid dynamics is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating fluid dynamics.");

DeclException1(
  MarangoniWithoutHeatTransferError,
  bool,
  << "Marangoni effect is activated (" << arg1
  << "), while heat transfer is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating heat transfer.");

DeclException1(
  MarangoniWithoutSurfaceTensionForceError,
  bool,
  << "Marangoni effect is activated (" << arg1
  << "), while surface tension force is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating surface tension force.");

DeclException1(
  SurfaceTensionForceWithoutVOFError,
  bool,
  << "Surface tension force is activated (" << arg1
  << "), while VOF is not activated (false)." << std::endl
  << "Surface tension force cannot be activated without activating VOF.");

DeclException1(
  MarangoniWithoutVOFError,
  bool,
  << "Marangoni effect is activated (" << arg1
  << "), while VOF is not activated (false)." << std::endl
  << "Marangoni effect cannot be activated without activating VOF.");

DeclException1(
  InterfaceSharpeningWithoutVOFError,
  bool,
  << "Interface sharpening is activated (" << arg1
  << "), while VOF is not activated (false)." << std::endl
  << "Interface sharpening cannot be activated without activating VOF.");

template <int dim>
MultiphysicsInterface<dim>::MultiphysicsInterface(
  const SimulationParameters<dim> &                            nsparam,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control,
  ConditionalOStream &               p_pcout)
  : multiphysics_parameters(nsparam.multiphysics)
  , verbosity(nsparam.non_linear_solver.verbosity)
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
      active_physics.push_back(PhysicsID::heat_transfer);
      physics[PhysicsID::heat_transfer] = std::make_shared<HeatTransfer<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }
  if (multiphysics_parameters.tracer)
    {
      active_physics.push_back(PhysicsID::tracer);
      physics[PhysicsID::tracer] = std::make_shared<Tracer<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }
  if (multiphysics_parameters.VOF)
    {
      active_physics.push_back(PhysicsID::VOF);
      physics[PhysicsID::VOF] = std::make_shared<VolumeOfFluid<dim>>(
        this, nsparam, p_triangulation, p_simulation_control);
    }
}

template <int dim>
TrilinosWrappers::MPI::Vector *
MultiphysicsInterface<dim>::get_projected_phase_fraction_gradient_solution()
{
  // Throw error if VOF is not enabled
  AssertThrow((std::find(active_physics.begin(),
                         active_physics.end(),
                         PhysicsID::VOF) != active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.vof_parameters.surface_tension_force.enable),
    ExcInternalError());

  return std::dynamic_pointer_cast<VolumeOfFluid<dim>>(physics[PhysicsID::VOF])
    ->get_projected_phase_fraction_gradient_solution();
}

template <int dim>
TrilinosWrappers::MPI::Vector *
MultiphysicsInterface<dim>::get_curvature_solution()
{
  // Throw error if VOF is not enabled
  AssertThrow((std::find(active_physics.begin(),
                         active_physics.end(),
                         PhysicsID::VOF) != active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.vof_parameters.surface_tension_force.enable),
    ExcInternalError());

  return std::dynamic_pointer_cast<VolumeOfFluid<dim>>(physics[PhysicsID::VOF])
    ->get_curvature_solution();
}

template <int dim>
DoFHandler<dim> *
MultiphysicsInterface<dim>::get_curvature_dof_handler()
{
  // Throw error if VOF is not enabled
  AssertThrow((std::find(active_physics.begin(),
                         active_physics.end(),
                         PhysicsID::VOF) != active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.vof_parameters.surface_tension_force.enable),
    ExcInternalError());

  return std::dynamic_pointer_cast<VolumeOfFluid<dim>>(physics[PhysicsID::VOF])
    ->get_curvature_dof_handler();
}

template <int dim>
DoFHandler<dim> *
MultiphysicsInterface<dim>::get_projected_phase_fraction_gradient_dof_handler()
{
  // Throw error if VOF is not enabled
  AssertThrow((std::find(active_physics.begin(),
                         active_physics.end(),
                         PhysicsID::VOF) != active_physics.end()),
              ExcInternalError());
  // Throw error if surface tension force is not enabled
  AssertThrow(
    (multiphysics_parameters.vof_parameters.surface_tension_force.enable),
    ExcInternalError());

  return std::dynamic_pointer_cast<VolumeOfFluid<dim>>(physics[PhysicsID::VOF])
    ->get_projected_phase_fraction_gradient_dof_handler();
}

template <int dim>
void
MultiphysicsInterface<dim>::inspect_multiphysics_models_dependencies(
  const SimulationParameters<dim> &nsparam)
{
  bool buoyancy_force_enabled = nsparam.multiphysics.buoyancy_force;
  bool heat_transfer_enabled  = nsparam.multiphysics.heat_transfer;
  bool marangoni_effect_enabled =
    nsparam.multiphysics.vof_parameters.surface_tension_force
      .enable_marangoni_effect;
  bool surface_tension_force_enabled =
    nsparam.multiphysics.vof_parameters.surface_tension_force.enable;
  bool fluid_dynamics_enabled = nsparam.multiphysics.fluid_dynamics;
  bool interface_sharpening_enabled =
    nsparam.multiphysics.vof_parameters.sharpening.enable;
  bool VOF_enabled = nsparam.multiphysics.VOF;

  // To avoid getting unused parameter warning
  _unused(buoyancy_force_enabled && heat_transfer_enabled &&
          marangoni_effect_enabled && surface_tension_force_enabled &&
          fluid_dynamics_enabled && interface_sharpening_enabled &&
          VOF_enabled);

  // Dependence of buoyant force on fluid dynamics
  Assert(!(buoyancy_force_enabled == true && fluid_dynamics_enabled == false),
         BuoyancyWithoutFluidDynamicsError(buoyancy_force_enabled));

  // Dependence of buoyant force on heat transfer
  Assert(!(buoyancy_force_enabled == true && heat_transfer_enabled == false),
         BuoyancyWithoutHeatTransferError(buoyancy_force_enabled));

  // Dependence of Marangoni effect on fluid dynamics
  Assert(!(marangoni_effect_enabled == true && fluid_dynamics_enabled == false),
         MarangoniWithoutFluidDynamicsError(marangoni_effect_enabled));

  // Dependence of Marangoni effect on heat transfer
  Assert(!(marangoni_effect_enabled == true && heat_transfer_enabled == false),
         MarangoniWithoutHeatTransferError(marangoni_effect_enabled));

  // Dependence of Marangoni effect on surface tension force
  Assert(!(marangoni_effect_enabled == true &&
           surface_tension_force_enabled == false),
         MarangoniWithoutSurfaceTensionForceError(marangoni_effect_enabled));

  // Dependence of surface tension force on VOF
  Assert(!(surface_tension_force_enabled == true && VOF_enabled == false),
         SurfaceTensionForceWithoutVOFError(surface_tension_force_enabled));

  // Dependence of Marangoni effect on VOF
  Assert(!(marangoni_effect_enabled == true && VOF_enabled == false),
         MarangoniWithoutVOFError(marangoni_effect_enabled));

  // Dependence of interface sharpening on VOF
  Assert(!(interface_sharpening_enabled == true && VOF_enabled == false),
         InterfaceSharpeningWithoutVOFError(interface_sharpening_enabled));
}

template class MultiphysicsInterface<2>;
template class MultiphysicsInterface<3>;
