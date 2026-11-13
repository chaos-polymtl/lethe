// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cfd_dem_simulation_parameters_h
#define lethe_cfd_dem_simulation_parameters_h

#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>

#include <solvers/simulation_parameters.h>

#include <dem/dem_solver_parameters.h>

#include <deal.II/base/exceptions.h>

#include <filesystem>
#include <string>

template <int dim>
class CFDDEMSimulationParameters
{
public:
  SimulationParameters<dim> cfd_parameters;
  DEMSolverParameters<dim>  dem_parameters;

  std::shared_ptr<Parameters::VoidFractionParameters<dim>> void_fraction;
  Parameters::CFDDEM                                       cfd_dem;

  void
  declare(ParameterHandler             &prm,
          Parameters::SizeOfSubsections size_of_subsections)
  {
    cfd_parameters.declare(prm, size_of_subsections);
    dem_parameters.declare(prm, size_of_subsections);

    void_fraction = std::make_shared<Parameters::VoidFractionParameters<dim>>();
    void_fraction->declare_parameters(prm);
    Parameters::CFDDEM::declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    cfd_parameters.parse(prm);
    dem_parameters.parse(prm);

    void_fraction->parse_parameters(prm);
    cfd_dem.parse_parameters(prm);

    // The checkpoints of the simulation must not overwrite the DEM checkpoint
    // it starts from.
    if (void_fraction->read_dem && cfd_parameters.restart_parameters.checkpoint)
      {
        const auto normalized_path = [](const std::string &path) {
          return std::filesystem::absolute(path).lexically_normal();
        };
        AssertThrow(
          normalized_path(cfd_parameters.simulation_control.output_folder +
                          cfd_parameters.restart_parameters.filename) !=
            normalized_path(void_fraction->dem_file_name),
          ExcMessage(
            "The checkpoints of this simulation, given by 'subsection "
            "simulation control' - 'set output path' and 'subsection restart' "
            "- 'set filename', would overwrite the DEM checkpoint it starts "
            "from, given by 'subsection void fraction' - 'set dem file name' (" +
            void_fraction->dem_file_name + "). Give them different names."));
      }
  }
};

#endif
