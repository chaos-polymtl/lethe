#ifndef lethe_euler_euler_sf_h
#define lethe_euler_euler_sf_h

#include <solvers/euler_euler_prm.h>
#include <solvers/euler_void_fraction.h>
#include <solvers/solid_phase.h>

#include <fem-dem/fluid_dynamics_vans.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <memory>

using namespace dealii;



template <int dim>
class EulerEulerOneWay : public FluidDynamicsVANS<dim>
{
public:
  EulerEulerOneWay(CFDDEMSimulationParameters<dim> &fluid_parameters,
                   const SolidPhaseParameters<dim> &solid_parameters);

  void
  solve() override;

private:
  SolidPhaseSolver<dim>       solid_solver;
  EulerEulerVoidFraction<dim> euler_void_fraction;

  void
  pass_fluid_solution_to_solid();

  void
  update_void_fraction_from_solid();

  GlobalVectorType fluid_drag_rhs;

  void
  assemble_fluid_drag_exchange_rhs();

  void
  inspect_mesh_boundaries() const;
};

#endif
