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

#include <core/boundary_conditions.h>
#include <core/rheological_model.h>
#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>

#ifndef lethe_isothermal_compressible_navier_stokes_assembler_h
#  define lethe_isothermal_compressible_navier_stokes_assembler_h


/**
 * @brief Class that assembles the core of the isothermal compressible Navier-Stokes equations.
 * This class assembles the weak form of:
 * $$\nabla \cdot (\rho \mathbf{u}) + \rho \mathbf{u} \cdot \nabla \mathbf{u} -
 * \nabla p - \mu \nabla^2 \mathbf{u} = 0 $$ with a full GLS stabilization
 * including the laplacian of the test function.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSIsothermalCompressibleNavierStokesAssemblerCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSIsothermalCompressibleNavierStokesAssemblerCore(
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

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the isothermal compressible Navier Stokes equations. For
 * example, if a BDF1 scheme is chosen, the following is assembled
 * $$\frac{(\rho \mathbf{u})^{t+\Delta t}-(\rho \mathbf{u})^{t}{\Delta t} +
 * \frac{(\psi \mathbf{p})^{t+\Delta t}-(\psi \mathbf{p})^{t}{\Delta t} $$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSIsothermalCompressibleNavierStokesAssemblerBDF
  : public NavierStokesAssemblerBase<dim>
{
public:
  GLSIsothermalCompressibleNavierStokesAssemblerBDF(
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
