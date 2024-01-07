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
 *
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
 */

/**
 * @brief This code tests the reynolds stress calculations in 2d with
 * Trilinos block vectors and MPI rank of 1 and 2 with adaptive time.
 */

// Deal.II includes
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/postprocessing_velocities.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // SimulationControl parameters
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.dt                           = 0.1;
  simulation_control_parameters.timeEnd                      = 1.0;
  simulation_control_parameters.adapt                        = true;
  simulation_control_parameters.adaptative_time_step_scaling = 0.99;

  Parameters::PostProcessing postprocessing_parameters;
  postprocessing_parameters.calculate_average_velocities = true;
  postprocessing_parameters.initial_time                 = 0.5;

  auto simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  std::vector<IndexSet> locally_owned_dofs(2);
  std::vector<IndexSet> locally_relevant_dofs(2);

  locally_owned_dofs[0].set_size(12);
  locally_owned_dofs[1].set_size(12);
  locally_owned_dofs[0].add_range(0, 8);
  locally_owned_dofs[1].add_range(0, 4);

  // Make triangulation and dummy dof_handler to construct average velocities
  parallel::distributed::Triangulation<2> tria(
    mpi_communicator,
    typename Triangulation<2>::MeshSmoothing(
      Triangulation<2>::smoothing_on_refinement |
      Triangulation<2>::smoothing_on_coarsening));
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<2> dof_handler(tria);


  AverageVelocities<2, GlobalBlockVectorType, std::vector<IndexSet>>
    postprocessing_velocities(dof_handler);

  GlobalBlockVectorType solution(locally_owned_dofs, mpi_communicator);
  solution.block(0)[0] = 2.0;
  solution.block(0)[1] = 0.1;
  solution.block(0)[2] = 2.5;
  solution.block(0)[3] = 0.56;
  solution.block(0)[4] = 0;
  solution.block(0)[5] = 7.9;
  solution.block(0)[6] = 1;
  solution.block(0)[7] = 1.2;
  solution.block(1)[0] = 30;
  solution.block(1)[1] = 20;
  solution.block(1)[2] = 26;
  solution.block(1)[2] = 15;

  GlobalBlockVectorType reynolds_normal_stresses;
  GlobalBlockVectorType reynolds_shear_stresses;

  // Time info
  const double time_end     = simulation_control_parameters.timeEnd;
  const double initial_time = postprocessing_parameters.initial_time;
  double       time         = simulation_control->get_current_time();
  double       dt           = 0.0;
  double       epsilon      = 1e-6;

  // Initialize averaged vectors
  postprocessing_velocities.initialize_vectors(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               2,
                                               mpi_communicator);

  // Time loop
  while (time < (time_end + epsilon)) // Until time reached end time
    {
      if (time > (initial_time - epsilon)) // Time reached the initial time
        {
          postprocessing_velocities.calculate_average_velocities(
            solution,
            postprocessing_parameters,
            simulation_control->get_current_time(),
            simulation_control->get_time_step());

          reynolds_normal_stresses =
            postprocessing_velocities.get_reynolds_normal_stresses();
          reynolds_shear_stresses =
            postprocessing_velocities.get_reynolds_shear_stresses();

          deallog << " Time  :      " << time << std::endl;
          deallog << " Time step  : " << dt << std::endl;
          deallog << " <u'u'> :     " << reynolds_normal_stresses.block(0)[0]
                  << " " << reynolds_normal_stresses.block(0)[2] << " "
                  << reynolds_normal_stresses.block(0)[4] << " "
                  << reynolds_normal_stresses.block(0)[6] << std::endl;
          deallog << " <v'v'> :     " << reynolds_normal_stresses.block(0)[1]
                  << " " << reynolds_normal_stresses.block(0)[3] << " "
                  << reynolds_normal_stresses.block(0)[5] << " "
                  << reynolds_normal_stresses.block(0)[7] << std::endl;
          deallog << " <u'v'> :     " << reynolds_shear_stresses.block(0)[0]
                  << " " << reynolds_shear_stresses.block(0)[2] << " "
                  << reynolds_shear_stresses.block(0)[4] << " "
                  << reynolds_shear_stresses.block(0)[6] << std::endl;
          deallog << " k :           " << reynolds_normal_stresses.block(1)[0]
                  << " " << reynolds_normal_stresses.block(1)[1] << " "
                  << reynolds_normal_stresses.block(1)[2] << " "
                  << reynolds_normal_stresses.block(1)[3] << std::endl;
          deallog << "" << std::endl;
        }

      // New solution values for next step
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
