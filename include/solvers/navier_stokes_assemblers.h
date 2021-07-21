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
#include <solvers/navier_stokes_scratch_data.h>

#ifndef lethe_navier_stokes_assemblers_h
#  define lethe_navier_stokes_assemblers_h

template <int dim>
class NavierStokesAssemblerBase
{
public:
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;

  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) = 0;
};



template <int dim>
class GLSNavierStokesAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerCore(
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
class GLSNavierStokesAssemblerSRF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerSRF(Parameters::VelocitySource velocity_sources)
    : velocity_sources(velocity_sources)
  {}

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  Parameters::VelocitySource velocity_sources;
};


template <int dim>
class GLSNavierStokesAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

template <int dim>
class GLSNavierStokesAssemblerSDIRK : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesAssemblerSDIRK(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

#endif
