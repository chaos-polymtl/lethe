// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
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
  /* This test checks the InterfaceTools::SignedDistanceSolver class from within
  an arbitrary structure that replicates the usual physical solver structure,
  denoted by the background prefix. This test performs verification of the
  signed distance computations for a sphere.
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
  for (int n = 0; n < dim; n++)
    {
      p_0[n] = 0.0;
      p_1[n] = 1.0;
    }
  GridGenerator::hyper_rectangle(*background_triangulation, p_0, p_1);
  background_triangulation->refine_global(3);

  // Discretize the domain as in an arbitrary background solver
  std::shared_ptr<FiniteElement<dim>> background_fe =
    std::make_shared<FE_Q<dim>>(1);
  std::shared_ptr<MappingQ<dim>> background_mapping =
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
  for (int n = 0; n < dim; n++)
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
  double max_reinitialization_distance = 1.0;
  std::shared_ptr<InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>
    signed_distance_solver = std::make_shared<
      InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>(
      background_triangulation,
      background_fe,
      max_reinitialization_distance,
      0.0,
      1.0,
      Parameters::Verbosity::quiet);

  signed_distance_solver->setup_dofs();
  signed_distance_solver->set_level_set_from_background_mesh(
    background_dof_handler, background_level_set);

  // Solve the signed_distance field.
  signed_distance_solver->solve();

  // Get the signed_distance field from the SignedDistanceSolver
  auto &signed_distance = signed_distance_solver->get_signed_distance();

  // Interpolate the signed_distance solution from the SignedDistanceSolver
  // DoFHandler to the main solver DoFHandler.
  GlobalVectorType background_signed_distance(background_locally_owned_dofs,
                                              background_locally_relevant_dofs,
                                              mpi_communicator);

  TrilinosWrappers::MPI::Vector tmp_background_signed_distance(
    background_locally_owned_dofs, mpi_communicator);

  FETools::interpolate(signed_distance_solver->dof_handler,
                       signed_distance,
                       background_dof_handler,
                       tmp_background_signed_distance);

  background_signed_distance = tmp_background_signed_distance;

  // Compute the L2 norm of the error.
  Vector<float> error_per_cell(background_triangulation->n_active_cells());

  VectorTools::integrate_difference(
    *background_mapping,
    background_dof_handler,
    background_signed_distance,
    Functions::SignedDistance::Sphere<dim>(sphere_center, sphere_radius),
    error_per_cell,
    QGauss<dim>(background_fe->degree + 1),
    VectorTools::L2_norm);

  const double error_L2 =
    VectorTools::compute_global_error(*background_triangulation,
                                      error_per_cell,
                                      VectorTools::L2_norm);

  deallog << "The L2 norm of the signed distance error in " << dim
          << "D is: " << error_L2 << std::endl;
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
