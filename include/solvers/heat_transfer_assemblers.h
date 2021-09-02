/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

 *
 * Implementation of heat transfer as an auxiliary physics.
 * This heat equation is weakly coupled to the velocity field.
 * Equation solved:
 * rho * Cp * (dT/dt + u.gradT) = k div(gradT) + nu/rho * (gradu : gradu)
 *
 * Author: Bruno Blais, Jeanne Joachim and Shahab Golshan, Polytechnique
 Montreal, 2020-
 */

#ifndef lethe_heat_transfer_assemblers_h
#define lethe_heat_transfer_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/heat_transfer_scratch_data.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for heat transfer
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerBase
{
public:
  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the heat transfer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &   copy_data) = 0;

  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the heat transfer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &   copy_data) = 0;
};

/**
 * @brief Class that assembles the core of the heat transfer equation.
 * This class assembles the weak form of:
 * $$ - k * \nabla^2 T + \rho * cp * \mathbf{u} * \nabla T - f - \tau :
 * \nabla \mathbf{u} =0 $$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class HeatTransferAssemblerCore : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::PhysicalProperties     physical_properties)
    : simulation_control(simulation_control)
    , physical_properties(physical_properties)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &   copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &   copy_data) override;

  const bool GGLS = true;

  std::shared_ptr<SimulationControl> simulation_control;
  Parameters::PhysicalProperties     physical_properties;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the heat transfer solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerBDF : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &   copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &   copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

/**
 * @brief Class that assembles Robin boundary condition, loop on faces (Newton's cooling law)
 * implementation similar to deal.ii step-7 for the heat transfer solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class HeatTransferAssemblerRBC : public HeatTransferAssemblerBase<dim>
{
public:
  HeatTransferAssemblerRBC(
    std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(HeatTransferScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &   copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(HeatTransferScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &   copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

#endif
