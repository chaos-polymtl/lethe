/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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
 */

#ifndef lethe_cfd_dem_simulation_parameters_h
#define lethe_cfd_dem_simulation_parameters_h

#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>

#include <solvers/analytical_solutions.h>
#include <solvers/initial_conditions.h>
#include <solvers/simulation_parameters.h>
#include <solvers/source_terms.h>

#include <dem/dem_solver_parameters.h>
#include <fem-dem/parameters_cfd_dem.h>

template <int dim>
class CFDDEMSimulationParameters
{
public:
  SimulationParameters<dim> cfd_parameters;
  DEMSolverParameters<dim>  dem_parameters;

  std::shared_ptr<Parameters::VoidFraction<dim>> void_fraction;
  Parameters::CFDDEM                             cfd_dem;

  void
  declare(ParameterHandler             &prm,
          Parameters::SizeOfSubsections size_of_subsections)
  {
    cfd_parameters.declare(prm, size_of_subsections);
    dem_parameters.declare(prm);

    void_fraction = std::make_shared<Parameters::VoidFraction<dim>>();
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
  }
};

#endif
