// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_assemblers_h
#define lethe_vof_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/physics_assemblers.h>
#include <solvers/vof_scratch_data.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the VOF solver.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
using VOFAssemblerBase =
  PhysicsAssemblerBase<VOFScratchData<dim>, StabilizedMethodsCopyData>;


/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for VOF subequation solvers.
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
using VOFSubequationAssemblerBase =
  PhysicsAssemblerBase<ScratchDataType, StabilizedMethodsCopyData>;

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
  VOFAssemblerCore(const std::shared_ptr<SimulationControl> &simulation_control,
                   const Parameters::FEM                     fem_parameters,
                   const Parameters::VOF                     vof_parameters)
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
  assemble_matrix(const VOFScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const VOFScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::FEM                    fem_parameters;
  const Parameters::VOF                    vof_parameters;

  // Controls if the compressibility term is assembled in the VOF equation
  const bool compressible;
};


/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the VOF equation. For example, if a BDF1 scheme is
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
  VOFAssemblerBDF(const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const VOFScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const VOFScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Assemble the Discontinuity-Capturing Directional Dissipation (DCDD)
 * stabilization term for the VOF phase fraction.
 *
 * @note For more information see Tezduyar, T. E. (2003). Computation of
 * moving boundaries and interfaces and stabilization parameters. International
 * Journal for Numerical Methods in Fluids, 43(5), 555-575. The implementation
 * is based on equations (70) and (79), which are adapted for the VOF solver.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAssemblerDCDDStabilization : public VOFAssemblerBase<dim>
{
public:
  /**
   * @brief Default constructor of the assembler.
   *
   * @param[in] simulation_control SimulationControl object that holds
   * information related to the control of the steady-state or transient
   * simulation. This is used to extrapolate velocity solutions in time for
   * transient simulations.
   */
  VOFAssemblerDCDDStabilization(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Assemble the matrix
   *
   * @param[in] scratch_data (see base class).
   *
   * @param[in,out] copy_data (see base class).
   */
  virtual void
  assemble_matrix(const VOFScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief Assemble the right-hand side (rhs).
   *
   * @param[in] scratch_data (see base class)
   *
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const VOFScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief VOF phase fraction gradient L2 projection assemblers.
 *
 * @tparam dim Number of dimensions of the problem.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAssemblerPhaseGradientProjection
  : public VOFSubequationAssemblerBase<
      VOFPhaseGradientProjectionScratchData<dim>>
{
public:
  /**
   * @brief Constructor of the assembler for phase fraction gradient (pfg)
   * L2 projection.
   *
   * @param[in] vof_parameters VOF simulation parameters.
   */
  VOFAssemblerPhaseGradientProjection(const Parameters::VOF &vof_parameters)
    : VOFSubequationAssemblerBase<VOFPhaseGradientProjectionScratchData<dim>>()
    , vof_parameters(vof_parameters)
  {}

  /**
   * @brief Default destructor.
   */
  ~VOFAssemblerPhaseGradientProjection() = default;

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data Scratch data containing the VOF phase gradient
   * projection information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   *
   * @param[in,out] copy_data Destination where the local_matrix is copied to.
   */
  void
  assemble_matrix(
    const VOFPhaseGradientProjectionScratchData<dim> &scratch_data,
    StabilizedMethodsCopyData                        &copy_data) override;

  /**
   * @brief Assemble the right-hand side (rhs).
   *
   * @param[in] scratch_data Scratch data containing the VOF phase gradient
   * projection information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   *
   * @param[in,out] copy_data Destination where the local_rhs is copied to.
   */
  void
  assemble_rhs(const VOFPhaseGradientProjectionScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

private:
  const Parameters::VOF vof_parameters;
};


/**
 * @brief VOF curvature L2 projection assemblers.
 *
 * @tparam dim Number of dimensions of the problem.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 *
 * @ingroup assemblers
 */
template <int dim>
class VOFAssemblerCurvatureProjection
  : public VOFSubequationAssemblerBase<VOFCurvatureProjectionScratchData<dim>>
{
public:
  /**
   * @brief Constructor of the assembler for curvature L2 projection.
   *
   * @param[in] vof_parameters VOF simulation parameters.
   */
  VOFAssemblerCurvatureProjection(const Parameters::VOF &vof_parameters)
    : VOFSubequationAssemblerBase<VOFCurvatureProjectionScratchData<dim>>()
    , vof_parameters(vof_parameters)
  {}

  /**
   * @brief Default destructor.
   */
  ~VOFAssemblerCurvatureProjection() = default;

  /**
   * @brief Assemble the matrix.
   *
   * @param[in] scratch_data Scratch data containing the VOF curvature
   * projection information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   *
   * @param[in,out] copy_data Destination where the local_matrix is copied to.
   */
  void
  assemble_matrix(const VOFCurvatureProjectionScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData &copy_data) override;

  /**
   * @brief Assemble the right-hand side (rhs).
   *
   * @param[in] scratch_data Scratch data containing the VOF curvature
   * projection information.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   *
   * @param[in,out] copy_data Destination where the local_rhs is copied to.
   */
  void
  assemble_rhs(const VOFCurvatureProjectionScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData &copy_data) override;

private:
  const Parameters::VOF vof_parameters;
};


/**
 * @brief TODO AMISHGA
 */
template <typename ScratchDataType>
class VOFSubequationAssemblerFactory
{
public:
  /**
   * @ brief TODO AMISHGA
   * @param[in] vof_parameters
   * @return
   */
  static std::shared_ptr<VOFSubequationAssemblerBase<ScratchDataType>>
  initialize_assembler(const Parameters::VOF & /*vof_parameters*/)
  {
    throw std::runtime_error("Unsupported ScratchDataType");
  }
};


/**
 * @brief TODO AMISHGA
 */
template <int dim>
class VOFSubequationAssemblerFactory<VOFPhaseGradientProjectionScratchData<dim>>
{
public:
  /**
   * @ brief TODO AMISHGA
   * @param[in] vof_parameters
   * @return
   */
  static std::shared_ptr<
    VOFSubequationAssemblerBase<VOFPhaseGradientProjectionScratchData<dim>>>
  initialize_assembler(const Parameters::VOF &vof_parameters)
  {
    return std::make_shared<VOFAssemblerPhaseGradientProjection<dim>>(
      vof_parameters);
  }
};

/**
 * @brief TODO AMISHGA
 */
template <int dim>
class VOFSubequationAssemblerFactory<VOFCurvatureProjectionScratchData<dim>>
{
public:
  /**
   * @ brief TODO AMISHGA
   * @param[in] vof_parameters
   * @return
   */
  static std::shared_ptr<
    VOFSubequationAssemblerBase<VOFCurvatureProjectionScratchData<dim>>>
  initialize_assembler(const Parameters::VOF &vof_parameters)
  {
    return std::make_shared<VOFAssemblerCurvatureProjection<dim>>(
      vof_parameters);
  }
};

#endif
