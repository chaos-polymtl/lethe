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



// template <int dim>
// struct EulerEulerGridParameters
// {
//   Point<dim> p1;
//   Point<dim> p2;

//   std::vector<unsigned int> subdivisions;

//   unsigned int global_refinement = 0;

//   std::string direction1 = "x";
//   std::string direction2 = "y";
// };

// template <int dim>
// void
// make_euler_euler_grid(parallel::distributed::Triangulation<dim> &tria,
//                       const EulerEulerGridParameters<dim> &grid_parameters);

// template <int dim>
// class EulerEulerOneWay
// {
// public:
//   EulerEulerOneWay(const EulerEulerMeshParameters<dim> &mesh_parameters,
//                    const SolidPhaseParameters          &solid_parameters,
//                    CFDDEMSimulationParameters<dim>     &fluid_parameters,
//                    const MPI_Comm                      &mpi_communicator,
//                    const bool                           verbose = true);

//   void
//   run();

// private:
//   void
//   solve_solid();

//   void
//   build_and_pass_alpha_f();

//   void
//   solve_fluid();

//   void
//   solve_solid_one_step();

//   void
//   solve_fluid_one_step();



//   const EulerEulerMeshParameters<dim> mesh_parameters;

//   MPI_Comm           mpi_communicator;
//   ConditionalOStream pcout;


//   const bool                 verbose;
//   const SolidPhaseParameters solid_parameters;

//   std::unique_ptr<SolidPhaseSolver<dim>> solid_solver;
//   FluidDynamicsVANS<dim>                 fluid_solver;
//   EulerEulerVoidFraction<dim>            void_fraction_solver;
// };


template <int dim>
class EulerEulerOneWay : public FluidDynamicsVANS<dim>
{
public:
  EulerEulerOneWay(CFDDEMSimulationParameters<dim> &fluid_parameters,
                   const SolidPhaseParameters      &solid_parameters);

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
};

#endif
