/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

/**
 * @brief This code tests averaging values in time with deal.II vectors.
 */

// Deal.II includes
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>


// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/postprocessing_velocities.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Parameters
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.dt      = 0.1;
  simulation_control_parameters.timeEnd = 1.0;
  simulation_control_parameters.adapt   = false;

  Parameters::PostProcessing postprocessing_parameters;
  postprocessing_parameters.calculate_average_velocities = true;
  postprocessing_parameters.initial_time                 = 0.5;

  auto simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  IndexSet locally_owned_dofs(8);
  IndexSet locally_relevant_dofs(8);

  locally_owned_dofs.add_range(0, 8);
  locally_relevant_dofs.add_range(0, 8);

  // Make triangulation and dummy dof_handler to construct average velocities
  parallel::distributed::Triangulation<3> tria(
    mpi_communicator,
    typename Triangulation<3>::MeshSmoothing(
      Triangulation<3>::smoothing_on_refinement |
      Triangulation<3>::smoothing_on_coarsening));
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<3> dof_handler(tria);


  AverageVelocities<3, LinearAlgebra::distributed::Vector<double>, IndexSet>
    average(dof_handler);

  LinearAlgebra::distributed::Vector<double> solution(locally_owned_dofs,
                                                      locally_relevant_dofs,
                                                      mpi_communicator);

  solution(0) = 0.0;
  solution(1) = 2.5;
  solution(2) = 10;
  solution(3) = 154.2;

  LinearAlgebra::distributed::Vector<double> average_solution;
  // average_solution.update_ghost_values();

  // Time and output info
  const double time_end     = simulation_control_parameters.timeEnd;
  const double initial_time = postprocessing_parameters.initial_time;
  double       time         = simulation_control->get_current_time();
  const double epsilon      = 1e-6;

  // Initialize averaged vectors
  average.initialize_vectors(locally_owned_dofs,
                             locally_relevant_dofs,
                             4,
                             mpi_communicator);

  // Time loop
  while (time < (time_end + epsilon)) // Until time reached end time
    {
      if (time > (initial_time - epsilon)) // Time reached the initial time
        {
          average.calculate_average_velocities(
            solution,
            postprocessing_parameters,
            simulation_control->get_current_time(),
            simulation_control->get_time_step());

          average_solution = average.get_average_velocities();

          deallog << " Time :             " << time << std::endl;
          deallog << " Average solution : " << average_solution[0] << " "
                  << average_solution[1] << " " << average_solution[2] << " "
                  << average_solution[3] << std::endl;
          deallog << "" << std::endl;
        }

      // new solution values for next step
      solution *= 0.9;

      // Integrate to get the next time
      simulation_control->integrate();

      // Break if the next time from integrate() is the same because
      // time will never get over the time end, but the average
      // velocities
      // at this time is wanted.
      if (abs(time - simulation_control->get_current_time()) < epsilon)
        break;

      time = simulation_control->get_current_time();
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test();
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
