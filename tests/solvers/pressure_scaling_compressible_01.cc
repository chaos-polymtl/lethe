// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the parsing of the simulation parameters rejects
 * a pressure scaling factor different from 1 when the density is not constant.
 * The isothermal compressible assemblers do not apply the pressure scaling
 * factor to the Jacobian, so this combination would produce a wrong Newton
 * correction. A pressure scaling factor with a constant density and a
 * non-constant density without pressure scaling must both remain accepted.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

// Lethe
#include <core/parameters.h>

#include <solvers/simulation_parameters.h>

#include <string>

// Tests
#include <../tests/tests.h>

using namespace dealii;

/**
 * @brief Parse a default set of simulation parameters in which only the
 * density model of the fluid and the pressure scaling factor are changed, and
 * report whether the parameters were accepted or rejected.
 *
 * @param[in] density_model Name of the density model of fluid 0.
 *
 * @param[in] pressure_scaling_factor Value of the pressure scaling factor of
 * the stabilization subsection.
 */
template <int dim>
void
check_parse(const std::string &density_model,
            const std::string &pressure_scaling_factor)
{
  // A fresh set of parameters is built for each check so that no state is
  // carried over from a previous parse.
  SimulationParameters<dim>     simulation_parameters;
  ParameterHandler              prm;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 1;

  simulation_parameters.declare(prm, size_of_subsections);

  prm.enter_subsection("physical properties");
  {
    prm.enter_subsection("fluid 0");
    {
      prm.set("density model", density_model);
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  prm.enter_subsection("stabilization");
  {
    prm.set("pressure scaling factor", pressure_scaling_factor);
  }
  prm.leave_subsection();

  deallog << "Density model: " << density_model
          << ", pressure scaling factor: " << pressure_scaling_factor << " -> ";
  try
    {
      simulation_parameters.parse(prm);
      deallog << "accepted" << std::endl;
    }
  catch (const std::exception &)
    {
      deallog << "rejected" << std::endl;
    }
}

template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  // Non-constant density with pressure scaling must be rejected
  check_parse<dim>("isothermal_ideal_gas", "10");

  // Non-constant density without pressure scaling must be accepted
  check_parse<dim>("isothermal_ideal_gas", "1");

  // Constant density with pressure scaling must be accepted
  check_parse<dim>("constant", "10");
}

int
main(int argc, char **argv)
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
