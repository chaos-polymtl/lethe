#include <solvers/multiphysics_interface.h>

#include "solvers/heat_transfer.h"


template <int dim>
MultiphysicsInterface<dim>::MultiphysicsInterface(
  const NavierStokesSolverParameters<dim> &                    nsparam,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  ConditionalOStream &                                         p_pcout)
  : multiphysics_parameters(nsparam.multiphysics)
{
  if (multiphysics_parameters.fluid_dynamics)
    {
      enabled_physics.push_back(PhysicsID::fluid_dynamics);
    }
  if (multiphysics_parameters.heat_transfer)
    {
      enabled_physics.push_back(PhysicsID::heat_transfer);
      physics[PhysicsID::heat_transfer] =
        std::make_shared(HeatTransfer<dim>(nsparam, p_triangulation, p_pcout));
    }
}
