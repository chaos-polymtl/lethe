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

#include <solvers/auxiliary_physics.h>
#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>

#ifndef lethe_isothermal_compressible_navier_stokes_vof_assembler_h
#  define lethe_isothermal_compressible_navier_stokes_vof_assembler_h


/**
 * @brief Class that assembles the core of the Navier-Stokes equation with
 * free surface using VOF modeling.
 * This class assembles the weak form of:
 * $$\nabla \cdot (\rho \mathbf{u}) + \rho \mathbf{u} \cdot \nabla
 * \mathbf{u} - \nabla p - \mu \nabla\cdot (\nabla \mathbf{u} +
 * (\nabla \mathbf{u})^T) = 0 $$ with a full GLS stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSIsothermalCompressibleNavierStokesVOFAssemblerCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSIsothermalCompressibleNavierStokesVOFAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    const SimulationParameters<dim> &  nsparam)
    : simulation_control(simulation_control)
    , vof_parameters(nsparam.multiphysics.vof_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::VOF              vof_parameters;
};

template <int dim>
class GLSIsothermalCompressibleNavierStokesVOFAssemblerBDF
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSIsothermalCompressibleNavierStokesVOFAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(NavierStokesScratchData<dim> &        scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(NavierStokesScratchData<dim> &        scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};


#endif
