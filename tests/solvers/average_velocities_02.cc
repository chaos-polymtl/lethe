// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests averaging values in time with Trilinos vectors with
 * adaptative time step scaling.
 */

// Deal.II includes
#include <deal.II/base/index_set.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/postprocessing_velocities.h>

// Tests
#include <../tests/tests.h>

template <int dim>
void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // SimulationControl parameters
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.dt = 0.1;
  simulation_control_parameters.end_control =
    Parameters::SimulationControl::EndControl::time;
  simulation_control_parameters.time_end                             = 1.0;
  simulation_control_parameters.output_iteration_frequency           = 1;
  simulation_control_parameters.time_step_adaptation_required        = true;
  simulation_control_parameters.adapt_with_cfl                       = true;
  simulation_control_parameters.adaptative_time_step_scaling         = 0.95;
  simulation_control_parameters.max_dt                               = 1e6;
  simulation_control_parameters.output_time_interval                 = {0,
                                                                        1.7976931348623157e3};
  simulation_control_parameters.time_step_independent_of_end_time    = true;
  simulation_control_parameters.adapt_with_capillary_time_step_ratio = false;

  Parameters::PostProcessing<dim> postprocessing_parameters;
  postprocessing_parameters.calculate_average_velocities        = true;
  postprocessing_parameters.initial_time_for_average_velocities = 0.5;

  auto simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  // Make triangulation and dof_handler to construct average velocities
  parallel::distributed::Triangulation<dim> tria(
    mpi_communicator,
    typename Triangulation<dim>::MeshSmoothing(
      Triangulation<dim>::smoothing_on_refinement |
      Triangulation<dim>::smoothing_on_coarsening));
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<dim> dof_handler(tria);

  // With a Q1 element on a single cell numbered component-wise, as in the
  // block solver, the velocity block holds (u, v, w) of the vertices, starting
  // with (0, 1, 2) for the first vertex, and the pressure block holds p of the
  // vertices, starting with 0 for the first vertex.
  dof_handler.distribute_dofs(FESystem<dim>(FE_Q<dim>(1), dim + 1));

  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  const std::vector<types::global_dof_index> dofs_per_block =
    DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
  const types::global_dof_index n_velocity_dofs = dofs_per_block[0];
  const types::global_dof_index n_dofs          = dof_handler.n_dofs();

  const IndexSet &owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet  relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  const std::vector<IndexSet> locally_owned_dofs = {
    owned_dofs.get_view(0, n_velocity_dofs),
    owned_dofs.get_view(n_velocity_dofs, n_dofs)};
  const std::vector<IndexSet> locally_relevant_dofs = {
    relevant_dofs.get_view(0, n_velocity_dofs),
    relevant_dofs.get_view(n_velocity_dofs, n_dofs)};

  AverageVelocities<dim, GlobalBlockVectorType, std::vector<IndexSet>> average(
    dof_handler);

  GlobalBlockVectorType solution(locally_owned_dofs, mpi_communicator);
  solution.block(0)[0] = 0.0;
  solution.block(0)[1] = 2.5;
  solution.block(0)[2] = 10;
  solution.block(1)[0] = 154.2;

  GlobalBlockVectorType average_solution(locally_owned_dofs, mpi_communicator);

  // Time and output info
  const double time_end = simulation_control_parameters.time_end;
  const double initial_time =
    postprocessing_parameters.initial_time_for_average_velocities;
  double       time    = simulation_control->get_current_time();
  double       dt      = 0.0;
  const double epsilon = 1e-6;

  // Initialize averaged vectors
  average.initialize_vectors(locally_owned_dofs,
                             locally_relevant_dofs,
                             mpi_communicator);

  // Time loop
  while (time < (time_end + dt)) // Until time reached end time
    {
      if (time > (initial_time - epsilon)) // Time reached the initial time
        {
          average.calculate_average_velocities(solution,
                                               postprocessing_parameters,
                                               time,
                                               dt);

          average_solution = *average.get_average_velocities();

          deallog << " Time :             " << time << std::endl;
          deallog << " Time step :        " << dt << std::endl;
          deallog << " Average solution : " << average_solution.block(0)[0]
                  << " " << average_solution.block(0)[1] << " "
                  << average_solution.block(0)[2] << " "
                  << average_solution.block(1)[0] << std::endl;
          deallog << "" << std::endl;
        }

      // new solution values for next step
      solution *= 0.9;

      // Integrate to get the next time
      simulation_control->integrate();

      // Break if the next time from integrate() is the same because
      // time will never get over the time end, but the average velocities
      // at this time is wanted.
      if (abs(time - simulation_control->get_current_time()) < epsilon)
        break;

      dt   = simulation_control->get_time_step();
      time = simulation_control->get_current_time();
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
