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
 */

#ifndef lethe_vof_assemblers_h
#define lethe_vof_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/vof_scratch_data.h>


/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the VOF solver
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAssemblerBase
{
public:
  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the VOF information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(VOFScratchData<dim>       &scratch_data,
                  StabilizedMethodsCopyData &copy_data) = 0;


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the VOF information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(VOFScratchData<dim>       &scratch_data,
               StabilizedMethodsCopyData &copy_data) = 0;
};


/**
 * @brief Class that assembles the core of the VOF solver.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot T \cdot \nabla T =0 \f$ with an SUPG
 * stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class VOFAssemblerCore : public VOFAssemblerBase<dim>
{
public:
  VOFAssemblerCore(std::shared_ptr<SimulationControl> simulation_control,
                   const Parameters::FEM              fem_parameters,
                   const Parameters::VOF              vof_parameters)
    : simulation_control(simulation_control)
    , fem_parameters(fem_parameters)
    , vof_parameters(vof_parameters)
    , compressible(vof_parameters.compressible)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(VOFScratchData<dim>       &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(VOFScratchData<dim>       &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::FEM              fem_parameters;
  const Parameters::VOF              vof_parameters;

  // Controls if the compressibility term is assembled in the VOF equation
  const bool compressible;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the VOF equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{T}^{t+\Delta t}-\mathbf{T}^{t}}{\Delta t}\f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAssemblerBDF : public VOFAssemblerBase<dim>
{
public:
  VOFAssemblerBDF(std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(VOFScratchData<dim>       &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(VOFScratchData<dim>       &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};

#endif
