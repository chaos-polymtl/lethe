// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_assemblers_h
#define lethe_cls_assemblers_h

#include <core/simulation_control.h>

#include <solvers/cls_scratch_data.h>
#include <solvers/copy_data.h>
#include <solvers/physics_assemblers.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the CLS solver.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
using CLSAssemblerBase =
  PhysicsAssemblerBase<CLSScratchData<dim>, StabilizedMethodsCopyData>;

/**
 * @brief A pure virtual class that serves as an interface for all
 * the assemblers for the CLS solver on internal faces.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
using CLSFaceAssemblerBase =
  PhysicsFaceAssemblerBase<CLSScratchData<dim>, StabilizedDGMethodsCopyData>;


/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for CLS subequation solvers.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 *
 *
 * @ingroup assemblers
 */
template <typename ScratchDataType>
using CLSSubequationAssemblerBase =
  PhysicsAssemblerBase<ScratchDataType, StabilizedMethodsCopyData>;

/**
 * @brief Class that assembles the core of the CLS solver.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot T \cdot \nabla T =0 \f$ with an SUPG
 * stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CLSAssemblerCore : public CLSAssemblerBase<dim>
{
public:
  CLSAssemblerCore(const std::shared_ptr<SimulationControl> &simulation_control,
                   const Parameters::FEM                     fem_parameters,
                   const Parameters::CLS                     cls_parameters)
    : simulation_control(simulation_control)
    , fem_parameters(fem_parameters)
    , cls_parameters(cls_parameters)
    , compressible(cls_parameters.compressible)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const CLSScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CLSScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::FEM                    fem_parameters;
  const Parameters::CLS                    cls_parameters;

  // Controls if the compressibility term is assembled in the CLS equation
  const bool compressible;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the CLS equation. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{T}^{t+\Delta t}-\mathbf{T}^{t}}{\Delta t}\f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CLSAssemblerBDF : public CLSAssemblerBase<dim>
{
public:
  CLSAssemblerBDF(const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const CLSScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CLSScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Assemble the Discontinuity-Capturing Directional Dissipation (DCDD)
 * stabilization term for the CLS phase indicator.
 *
 * @note For more information see Tezduyar, T. E. (2003). Computation of
 * moving boundaries and interfaces and stabilization parameters. International
 * Journal for Numerical Methods in Fluids, 43(5), 555-575. The implementation
 * is based on equations (70) and (79), which are adapted for the CLS solver.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class CLSAssemblerDCDDStabilization : public CLSAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor of the assembler.
   *
   * @param[in] simulation_control SimulationControl object that holds
   * information related to the control of the steady-state or transient
   * simulation. This is used to extrapolate velocity solutions in time for
   * transient simulations.
   *
   * @param[in] stabilization_parameters Parameters for the stabilization to
   * get the value of the diffusion coefficient.
   */
  CLSAssemblerDCDDStabilization(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::Stabilization          &stabilization_parameters)
    : simulation_control(simulation_control)
    , diffusion_constant(stabilization_parameters.dcdd_diffusion_coeff)
  {}

  /**
   * @brief Assemble the matrix
   *
   * @param[in] scratch_data (see base class).
   *
   * @param[in,out] copy_data (see base class).
   */
  virtual void
  assemble_matrix(const CLSScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief Assemble the right-hand side (rhs).
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CLSScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;

  /// Diffusion coefficient
  const double diffusion_constant;
};


/**
 * @brief Class that assembles the core (cells) of the CLS equation for DG elements.
 * This class assembles the weak form of:
 * \f$\mathbf{u} \cdot \nabla \phi = 0 \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CLSAssemblerDGCore : public CLSAssemblerBase<dim>
{
public:
  CLSAssemblerDGCore()
  {}

  /**
   * @brief Assemble the matrix
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const CLSScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief Assemble the rhs
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CLSScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;
};


/**
 * @brief Assembles the symmetric interior penalty Galerkin (SIPG) method (or
 * Nitsche's method) for internal faces. This assembler is only required
 * when solving the CLS equation using a discontinuous Galerkin
 * discretization.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CLSAssemblerSIPG : public CLSFaceAssemblerBase<dim>
{
public:
  CLSAssemblerSIPG()
  {}

  /**
   * @brief Interface to call the matrix assembly
   * @param[in]  scratch_data Scratch data containing the CLS
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied.
   */
  virtual void
  assemble_matrix(const CLSScratchData<dim>   &scratch_data,
                  StabilizedDGMethodsCopyData &copy_data) override;


  /**
   * @brief Interface for the call to rhs
   * @param[in]  scratch_data Scratch data containing the CLS
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied.
   */
  virtual void
  assemble_rhs(const CLSScratchData<dim>   &scratch_data,
               StabilizedDGMethodsCopyData &copy_data) override;
};


#endif
