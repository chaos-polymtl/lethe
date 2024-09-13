// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_assemblers_h
#define lethe_tracer_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/tracer_scratch_data.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Tracer equation
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerBase
{
public:
  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(TracerScratchData<dim>    &scratch_data,
                  StabilizedMethodsCopyData &copy_data) = 0;


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(TracerScratchData<dim>    &scratch_data,
               StabilizedMethodsCopyData &copy_data) = 0;
};


/**
 * @brief Class that assembles the core of the Tracer equation.
 * This class assembles the weak form of:
 * \f$\mathbf{u} \cdot \nabla T - D \nabla^2 =0 \f$ with an SUPG
 * stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class TracerAssemblerCore : public TracerAssemblerBase<dim>
{
public:
  TracerAssemblerCore(std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(TracerScratchData<dim>    &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(TracerScratchData<dim>    &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};



/**
 * @brief Class that assembles the core of the Tracer equation for DG elements.
 * This class assembles the weak form of:
 * \f$\mathbf{u} \cdot \nabla T - D \nabla^2 =0 \f$
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class TracerAssemblerDGCore : public TracerAssemblerBase<dim>
{
public:
  TracerAssemblerDGCore(std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(TracerScratchData<dim>    &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(TracerScratchData<dim>    &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Tracer equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{T}^{t+\Delta t}-\mathbf{T}^{t}}{\Delta t}\f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerBDF : public TracerAssemblerBase<dim>
{
public:
  TracerAssemblerBDF(std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(TracerScratchData<dim>    &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(TracerScratchData<dim>    &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief A pure virtual class that serves as an interface for boundary and internal faces that occurs when using a discontinuous Galerkin discretization.
 * The main difference between the TracerFaceAssembler and the TracerAssembler
 * is that the TracerFaceAssembler assembles the matrix and rhs for internal
 * faces and thus requires the StabilizedDGMethodsCopyData class.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerFaceAssembler
{
public:
  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(TracerScratchData<dim>      &scratch_data,
                  StabilizedDGMethodsCopyData &copy_data) = 0;


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(TracerScratchData<dim>      &scratch_data,
               StabilizedDGMethodsCopyData &copy_data) = 0;
};


/**
 * @brief Assembles the symmetric interior penalty (SIPG) method (or Nitche's method) for internal faces. This assembler is only required when solving the Tracer equation using a discontinuous Galerkin discretization.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerSIPG : public TracerFaceAssembler<dim>
{
public:
  TracerAssemblerSIPG()
  {}

  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_matrix(TracerScratchData<dim>      &scratch_data,
                  StabilizedDGMethodsCopyData &copy_data);


  /**
   * @brief assemble_matrix Interface for the call to rhs
   * @param scratch_data Scratch data containing the Tracer information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and loc
   */

  virtual void
  assemble_rhs(TracerScratchData<dim>      &scratch_data,
               StabilizedDGMethodsCopyData &copy_data);
};



#endif
