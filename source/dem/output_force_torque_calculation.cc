#include <dem/output_force_torque_calculation.h>

void
write_forces_torques_output_locally(
  std::map<unsigned int, Tensor<1, 3>> force_on_walls,
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls)
{
  TableHandler table;

  for (const auto &it : force_on_walls)
    {
      table.add_value("B_id", it.first);
      table.add_value("Fx", force_on_walls[it.first][0]);
      table.add_value("Fy", force_on_walls[it.first][1]);
      table.add_value("Fz", force_on_walls[it.first][2]);
      table.add_value("Tx", torque_on_walls[it.first][0]);
      table.add_value("Ty", torque_on_walls[it.first][1]);
      table.add_value("Tz", torque_on_walls[it.first][2]);
    }
  table.set_precision("Fx", 9);
  table.set_precision("Fy", 9);
  table.set_precision("Fz", 9);
  table.set_precision("Tx", 9);
  table.set_precision("Ty", 9);
  table.set_precision("Tz", 9);

  table.write_text(std::cout);
}

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
