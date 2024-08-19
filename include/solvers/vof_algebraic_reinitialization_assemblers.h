/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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

#ifndef lethe_vof_algebraic_reinitialization_assemblers_h
#define lethe_vof_algebraic_reinitialization_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/vof_algebraic_reinitialization_scratch_data.h>

/**
 * @brief Abstract class that serves as a base for all the assemblers of the
 * VOF Algebraic Reinitialization solver.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAlgebraicReinitializationAssemblerBase
{
public:
  /**
   * @brief Assembles matrix
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_matrix(
    const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
    StabilizedMethodsCopyData                          &copy_data) = 0;

  /**
   * @brief Assembles right-hand side (rhs)
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_rhs(const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) = 0;
};


/**
 * @brief Assembles the core of the algebraic interface reinitialization partial
 * differential equation.
 * The weak form of the core terms goes as follows:
 * \f$ \int_\Omega \nabla v \cdot \left[\epsilon (\nabla \phi_\mathrm{reinit}
 * \cdot \mathbf{n}) \mathbf{n} - (\phi_\mathrm{reinit} (1-\phi_\mathrm{reinit)
 * \mathbf{n}) \right] \mathrm{d} \Omega \f$ with \f$ v \f$ the test function.
 *
 * @remark Time integration term is handled by the
 * VOFAlgebraicReinitializationAssemblerBDF class.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAlgebraicReinitializationAssemblerCore
  : public VOFAlgebraicReinitializationAssemblerBase<dim>
{
public:
  VOFAlgebraicReinitializationAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    const Parameters::FEM              fem_parameters,
    const Parameters::VOF vof_algebraic_reinitialization_parameters)
    : simulation_control(simulation_control)
    , fem_parameters(fem_parameters)
    , vof_algebraic_reinitialization_parameters(
        vof_algebraic_reinitialization_parameters)
  {}

  /**
   * @brief Assembles matrix
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_matrix(
    const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
    StabilizedMethodsCopyData                          &copy_data) override;
  /**
   * @brief Assembles right-hand side (rhs)
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_rhs(const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::FEM              fem_parameters;
  const Parameters::VOF              vof_algebraic_reinitialization_parameters;
};


/**
 * @brief Assembles the time-stepping term of the algebraic interface
 * reinitialization partial differential equation (PDE).
 *
 * @note For example, if a BDF1 scheme is  chosen, the following is assembled:
 * \f$\int_\Omega v \frac{\phi_\mathrm{reinit}^{t+\Delta t}-
 * \phi_\mathrm{reinit}^{t}}{\Delta t} \mathrm{d}\Omega \f$ with \f$ v \f$ the
 * test function.
 *
 * @remark The other terms of the PDE is handled by the
 * VOFAlgebraicReinitializationAssemblerCore class.
 *
 * @remark At the moment, only BDF1 is implemented.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAlgebraicReinitializationAssemblerBDF
  : public VOFAlgebraicReinitializationAssemblerBase<dim>
{
public:
  VOFAlgebraicReinitializationAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Assembles matrix
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_matrix(
    const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
    StabilizedMethodsCopyData                          &copy_data) override;
  /**
   * @brief Assembles right-hand side (rhs)
   *
   * @param[in] scratch_data Scratch data containing the VOF algebraic interface
   * reinitialization information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix and right-hand side assembly.
   *
   * @param[in,out] copy_data Destination where the element for the local_rhs
   * and local_matrix are copied to.
   */
  virtual void
  assemble_rhs(const VOFAlgebraicReinitializationScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

#endif
