// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Check the output file of the force and torque exerted by the particles
 * on a solid object. The file of a new simulation contains a header followed by
 * one row per output time. When the simulation is restarted, the rows written
 * after the checkpoint are removed before new rows are appended. When the file
 * of a restarted simulation does not exist, it is created.
 */

// Deal.II
#include <deal.II/base/tensor.h>

// Lethe
#include <dem/output_force_torque_calculation.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <fstream>
#include <string>

using namespace dealii;

/**
 * @brief Print the content of a file in the log.
 *
 * @param[in] filename Name of the file.
 */
void
print_file(const std::string &filename)
{
  std::ifstream file(filename);
  std::string   line;
  while (std::getline(file, line))
    deallog << line << std::endl;
}

void
test()
{
  const std::string filename  = "solid_forces_00.dat";
  const double      time_step = 1e-3;

  // Tolerance used to find the rows written after the checkpoint
  const double time_tolerance = 0.5 * time_step;

  // New simulation, with three output times
  initialize_solid_forces_torques_file(filename, false, 0., time_tolerance);
  for (unsigned int i = 1; i <= 3; ++i)
    append_solid_forces_torques_to_file(
      filename,
      i * time_step,
      Tensor<1, 3>({1. * i, -2. * i, 0.5 * i}),
      Tensor<1, 3>({-0.1 * i, 0.2 * i, 0.3 * i}));

  deallog << "New simulation" << std::endl;
  print_file(filename);

  // Simulation restarted from a checkpoint written at the second output time.
  // The third row must be replaced by the one written after the restart.
  initialize_solid_forces_torques_file(filename,
                                       true,
                                       2 * time_step,
                                       time_tolerance);
  append_solid_forces_torques_to_file(filename,
                                      3 * time_step,
                                      Tensor<1, 3>({7., 8., 9.}),
                                      Tensor<1, 3>({-7., -8., -9.}));

  deallog << "Restarted simulation" << std::endl;
  print_file(filename);

  // Restarted simulation for which no file exists yet
  const std::string new_filename = "solid_forces_01.dat";
  initialize_solid_forces_torques_file(new_filename,
                                       true,
                                       2 * time_step,
                                       time_tolerance);
  append_solid_forces_torques_to_file(new_filename,
                                      3 * time_step,
                                      Tensor<1, 3>({1., 2., 3.}),
                                      Tensor<1, 3>({4., 5., 6.}));

  deallog << "Restarted simulation without an existing file" << std::endl;
  print_file(new_filename);

  // New simulation overwriting an existing file
  initialize_solid_forces_torques_file(filename, false, 0., time_tolerance);

  deallog << "New simulation overwriting an existing file" << std::endl;
  print_file(filename);
}

int
main()
{
  try
    {
      initlog();
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
