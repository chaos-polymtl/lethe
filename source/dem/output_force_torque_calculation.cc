// SPDX-FileCopyrightText: Copyright (c) 2021-2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/output_force_torque_calculation.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace
{
  /// Number of digits after the decimal point of the values written in the
  /// files of the force and torque on the solid objects. It is large enough to
  /// distinguish the time of consecutive iterations of long DEM simulations.
  constexpr unsigned int solid_forces_torques_precision = 12;

  /// Width of the columns of the files of the force and torque on the solid
  /// objects. It fits a signed value in scientific notation with the precision
  /// above, and leaves at least one blank between two values.
  constexpr unsigned int solid_forces_torques_column_width =
    solid_forces_torques_precision + 10;

  /// Names of the columns of the files of the force and torque on the solid
  /// objects.
  const std::array<std::string, 7> solid_forces_torques_column_names = {
    {"time", "f_x", "f_y", "f_z", "T_x", "T_y", "T_z"}};
} // namespace

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
  const std::string                               &filename,
  const unsigned int                               output_frequency,
  const std::vector<unsigned int>                 &boundary_index,
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

void
initialize_solid_forces_torques_file(const std::string &filename,
                                     const bool         restart,
                                     const double       restart_time,
                                     const double       time_tolerance)
{
  if (restart && std::filesystem::exists(filename))
    {
      // The rows are written in chronological order after the header. Find the
      // end of the last row written before the checkpoint, and discard what
      // follows. A last row without an end of line, which was interrupted
      // while being written, is also discarded since tellg() fails when the
      // end of the file is reached.
      std::ifstream  existing_file(filename);
      std::string    line;
      std::streamoff kept_size = 0;
      bool           is_header = true;
      while (std::getline(existing_file, line))
        {
          if (!is_header)
            {
              std::istringstream row(line);
              double             time;
              if (!(row >> time) || time > restart_time + time_tolerance)
                break;
            }

          const std::streampos end_of_line = existing_file.tellg();
          if (end_of_line == std::streampos(-1))
            break;

          kept_size = end_of_line;
          is_header = false;
        }
      existing_file.close();

      std::filesystem::resize_file(filename,
                                   static_cast<std::uintmax_t>(kept_size));

      // The existing file is kept if it has a header
      if (kept_size > 0)
        return;
    }

  std::ofstream file(filename);
  for (const auto &column_name : solid_forces_torques_column_names)
    file << std::setw(solid_forces_torques_column_width) << column_name;
  file << '\n';
}

void
append_solid_forces_torques_to_file(const std::string  &filename,
                                    const double        time,
                                    const Tensor<1, 3> &force,
                                    const Tensor<1, 3> &torque)
{
  std::ofstream file(filename, std::ios::app);
  file << std::scientific << std::setprecision(solid_forces_torques_precision);

  file << std::setw(solid_forces_torques_column_width) << time;
  for (unsigned int d = 0; d < 3; ++d)
    file << std::setw(solid_forces_torques_column_width) << force[d];
  for (unsigned int d = 0; d < 3; ++d)
    file << std::setw(solid_forces_torques_column_width) << torque[d];
  file << '\n';
}
