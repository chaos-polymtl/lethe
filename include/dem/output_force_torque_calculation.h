/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Lethe's community, 2021
 */

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/table_handler.h>

#ifndef lethe_output_force_torque_calculation_h
#  define lethe_output_force_torque_calculation_h

/**
 * @brief write_forces_torques_output_locally
 * Writes the results of force and torque calculations in the terminal at the
 * frequency requested in the prm.
 */
void
write_forces_torques_output_locally(
  std::map<unsigned int, Tensor<1, 3>> force_on_walls,
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls);

/**
 * @brief write_forces_torques_output_results
 * Writes the results of force and torque calculations in a file, and it depends
 * on the verbosity, in the terminal
 */
template <int dim>
void
write_forces_torques_output_results(
  const std::string                                filename,
  const unsigned int                               output_frequency,
  const std::vector<unsigned int>                  boundary_index,
  const double                                     time_step,
  DEM::dem_data_structures<3>::vector_on_boundary &forces_boundary_information,
  DEM::dem_data_structures<3>::vector_on_boundary &torques_boundary_information)
{
  unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  for (unsigned int i = 0; i < boundary_index.size(); i++)
    {
      std::string filename_boundary_description =
        std::to_string(boundary_index[i]);
      std::string filename_output =
        filename + "_boundary_" + filename_boundary_description;
      if (this_mpi_process == i || (i + 1) > n_mpi_processes)
        {
          TableHandler table;

          for (unsigned int j = 0; j < forces_boundary_information.size();
               j += output_frequency)
            {
              table.add_value("Force x", forces_boundary_information[j][i][0]);
              table.add_value("Force y", forces_boundary_information[j][i][1]);
              table.add_value("Force z", forces_boundary_information[j][i][2]);
              table.add_value("Torque x",
                              torques_boundary_information[j][i][0]);
              table.add_value("Torque y",
                              torques_boundary_information[j][i][1]);
              table.add_value("Torque z",
                              torques_boundary_information[j][i][2]);
              table.add_value("Current time", j * time_step);
            }

          table.set_precision("Force x", 9);
          table.set_precision("Force y", 9);
          table.set_precision("Force z", 9);
          table.set_precision("Torque x", 9);
          table.set_precision("Torque y", 9);
          table.set_precision("Torque z", 9);
          table.set_precision("Current time", 9);

          // output
          std::ofstream out_file(filename_output);
          table.write_text(out_file);
          out_file.close();
        }
    }
}

#endif // lethe_output_force_torque_calculation_h
