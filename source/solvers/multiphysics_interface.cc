#include <solvers/heat_transfer.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/tracer.h>
#include <solvers/vof.h>


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

template class MultiphysicsInterface<2>;
template class MultiphysicsInterface<3>;
