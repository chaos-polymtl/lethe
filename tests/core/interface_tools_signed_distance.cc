// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/interface_tools.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  /* This test checks the instanciation of the
  InterfaceTools::SignedDistanceSolver class from within an abitrary structure
  that replicates the usual physical solver structure, denoted by the background
  prefix. This is a dummy test intended to check the architecture, and not the
  SignedDistanceSolver performances yet! When the implementation of the
  SignedDistanceSolver will be completed, this test will be extended to a
  verification of the signed distance computation for a sphere.
  */
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Triangulation (as a shared_ptr to reproduce an arbitrary background solver
  // architecture)
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
    background_triangulation =
      std::make_shared<parallel::distributed::Triangulation<dim>>(
        mpi_communicator,
        typename Triangulation<dim>::MeshSmoothing(
          Triangulation<dim>::smoothing_on_refinement |
          Triangulation<dim>::smoothing_on_coarsening));

  Point<dim> p_0 = Point<dim>();
  Point<dim> p_1 = Point<dim>();
  for (unsigned int n = 0; n < dim; n++)
    {
      p_0[n] = 0.0;
      p_1[n] = 1.0;
    }
  GridGenerator::hyper_rectangle(*background_triangulation, p_0, p_1);
  background_triangulation->refine_global(3);

  // Discretize the domain as in an arbitrary background solver
  std::shared_ptr<FiniteElement<dim>> background_fe =
    std::make_shared<FE_Q<dim>>(1);
  std::shared_ptr<Mapping<dim>> background_mapping =
    std::make_shared<MappingQ<dim>>(background_fe->degree);

  DoFHandler<dim> background_dof_handler =
    DoFHandler(*background_triangulation);
  background_dof_handler.distribute_dofs(*background_fe);

  IndexSet background_locally_owned_dofs =
    background_dof_handler.locally_owned_dofs();
  IndexSet background_locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(background_dof_handler);

  // Set the background level-set field of the sphere
  GlobalVectorType background_level_set(background_locally_owned_dofs,
                                        background_locally_relevant_dofs,
                                        mpi_communicator);

  const double sphere_radius = 0.25;
  Point<dim>   sphere_center = Point<dim>();
  for (unsigned int n = 0; n < dim; n++)
    {
      sphere_center[n] = 0.5;
    }
  VectorTools::interpolate(
    *background_mapping,
    background_dof_handler,
    Functions::SignedDistance::Sphere<dim>(sphere_center, sphere_radius),
    background_level_set);

  // Instanciate the SignedDistanceSolver as it would as a member of the
  // background solver class
  double max_reinitialization_distance = 0.2;
  std::shared_ptr<InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>
    signed_distance_solver = std::make_shared<
      InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>(
      background_triangulation,
      background_fe,
      max_reinitialization_distance,
      0.0);

  signed_distance_solver->setup_dofs();
  signed_distance_solver->set_level_set_from_background_mesh(
    background_dof_handler, background_level_set);

  // Solve the signed_distance field. The solve() method only initialize the
  // signed_distance to the max_reinitialization_distance for now.
  signed_distance_solver->solve();

  // Get the signed_distance field from the SignedDistanceSolver
  auto &signed_distance = signed_distance_solver->get_signed_distance();

  // The extrema of signed_distance should be the +/-
  // max_reinitialization_distance. This is a dummy test intended to check the
  // architecture, and not the SignedDistanceSolver performances yet!
  deallog << "The signed distance interval in " << dim << "D is: ["
          << signed_distance.min() << "," << signed_distance.max() << "]"
          << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
