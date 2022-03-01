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


/**
 * @brief Class that assembles the core of the Navier-Stokes equation with
 * free surface using VOF modeling.
 * This class assembles the weak form of:
 * $$\rho \mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \mu \nabla^2
 * \mathbf{u} =0 $$ with an SUPG and PSPG stabilziation
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerCore : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerCore(
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

  const bool SUPG = true;

  std::shared_ptr<SimulationControl> simulation_control;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier-Stokes equation with
 * free surface using VOF modeling.. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * $$\frac{(\rho \mathbf{u})^{t+\Delta t}-(\rho \mathbf{u})^{t}{\Delta t}
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */

template <int dim>
class GLSNavierStokesVOFAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerBDF(
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
 * @brief Class that assembles the Continuum Surface Force (CSF) for the
 * Navier-Stokes equations. The following equation is assembled
 *
 * $$\mathbf{F_{CSV}}=\sigma \rho k \nabla \phi
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class GLSNavierStokesVOFAssemblerCSF : public NavierStokesAssemblerBase<dim>
{
public:
  GLSNavierStokesVOFAssemblerCSF(
    std::shared_ptr<SimulationControl> p_simulation_control,
    Parameters::SurfaceTensionForce    p_STF_properties)
    : simulation_control(p_simulation_control)
    , STF_properties(p_STF_properties)
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

  std::shared_ptr<SimulationControl>    simulation_control;
  const Parameters::SurfaceTensionForce STF_properties;
};

#endif
