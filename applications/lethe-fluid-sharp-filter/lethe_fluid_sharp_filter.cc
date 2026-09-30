// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/shape.h>
#include <core/utilities.h>

#include <fem-dem/anderson_jackson_filter.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/fluid_dynamics_sharp.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/timer.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace dealii;

namespace
{
  /**
   * @brief Restore a lethe-fluid-sharp snapshot, filter it with the
   * Anderson-Jackson filter and write the filtered fields.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] file_name Name of the parameter file.
   *
   * @param[in] pcout Output stream of the first process.
   */
  template <int dim>
  void
  run(const std::string &file_name, const ConditionalOStream &pcout)
  {
    ParameterHandler                  prm;
    CFDDEMSimulationParameters<dim>   NSparam;
    Parameters::AndersonJacksonFilter filter_parameters;
    NSparam.declare(prm, Parameters::get_size_of_subsections(file_name));
    Parameters::AndersonJacksonFilter::declare_parameters(prm);

    prm.parse_input(file_name);
    NSparam.parse(prm);
    filter_parameters.parse_parameters(prm);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        print_comment_to_output_file(pcout, prm);
        print_parameters_to_output_file(pcout, prm, file_name);
      }

    // The filtered fields must not overwrite the output of the simulation. The
    // output files are named by concatenating the folder and the name.
    const auto &simulation_control = NSparam.cfd_parameters.simulation_control;
    AssertThrow(
      std::filesystem::weakly_canonical(filter_parameters.output_folder +
                                        filter_parameters.output_name) !=
        std::filesystem::weakly_canonical(simulation_control.output_folder +
                                          simulation_control.output_name),
      ExcMessage("The output folder and name of the "
                 "Anderson-Jackson filter must differ from those "
                 "of the simulation."));

    FluidDynamicsSharp<dim>    problem(NSparam);
    AndersonJacksonFilter<dim> filter(filter_parameters,
                                      NSparam.cfd_parameters.timer,
                                      MPI_COMM_WORLD);
    {
      TimerOutput::Scope t(filter.get_computing_timer(), "Load restart");
      problem.initialize_for_postprocessing();
    }

    std::vector<std::shared_ptr<Shape<dim>>> shapes;
    for (const auto &particle : problem.get_particles())
      shapes.push_back(particle.shape);
    const ShapesImmersedSolidClassifier<dim> classifier(shapes);

    // Velocity imposed on the walls, at the time of the snapshot.
    const auto wall_velocities =
      make_wall_velocities(NSparam.cfd_parameters.boundary_conditions,
                           problem.get_current_time());

    filter.apply(problem.get_dof_handler(),
                 problem.get_fluid_mapping(),
                 problem.get_fluid_solution(),
                 classifier,
                 NSparam.cfd_parameters.boundary_conditions.periodic_boundaries,
                 wall_velocities);
    filter.write_output(problem.get_fluid_mapping(),
                        problem.get_current_time(),
                        simulation_control.group_files);
    filter.print_timer_summary();
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      ConditionalOStream pcout(
        std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

      auto [options, args] = parse_args(argc, argv);

      // Print version information
      if (options["-V"])
        {
          pcout << "Running: " << concatenate_strings(argc, argv) << std::endl;

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_version_info(pcout);

          return EXIT_SUCCESS;
        }

      if (args.empty())
        {
          pcout << "Usage: " << argv[0] << " input_file" << std::endl;
          return EXIT_FAILURE;
        }

      const std::string  file_name(args[0]);
      const unsigned int dim = get_dimension(file_name);

      if (dim == 2)
        run<2>(file_name, pcout);
      else if (dim == 3)
        run<3>(file_name, pcout);
      else
        return EXIT_FAILURE;
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
