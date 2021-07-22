/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * ---------------------------------------------------------------------*/


#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>

#ifndef lethe_navier_stokes_vof_assemblers_h
#  define lethe_navier_stokes_vof_assemblers_h


template <int dim>
class GLSNavierStokesFreeSurfaceAssemblerCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesFreeSurfaceAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
  {}


  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const bool SUPG = true;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
};

template <int dim>
class GLSNavierStokesFreeSurfaceAssemblerBDF
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesFreeSurfaceAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
  {}

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
};

#endif
