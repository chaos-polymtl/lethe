// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_assemblers_h
#define lethe_tracer_assemblers_h

#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/physics_assemblers.h>
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
using TracerAssemblerBase =
  PhysicsAssemblerBase<TracerScratchData<dim>, StabilizedMethodsCopyData>;

/**
 * @brief A pure virtual class that serves as an interface for boundary and
 * internal faces that occur when using a discontinuous Galerkin discretization.
 * The main difference between the TracerFaceAssembler and the TracerAssembler
 * is that the TracerFaceAssembler assembles the matrix and rhs for internal
 * faces and thus requires the StabilizedDGMethodsCopyData class.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
using TracerFaceAssemblerBase =
  PhysicsFaceAssemblerBase<TracerScratchData<dim>, StabilizedDGMethodsCopyData>;

/**
 * @brief Class that assembles the core of the Tracer equation.
 * This class assembles the weak form of:
 * \f$\mathbf{u} \cdot \nabla T - D \nabla^2T = S \f$ with an SUPG
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
  TracerAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Assembles the matrix
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;

  /**
   * @brief Assembles the rhs
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the core (cells) of the Tracer equation for DG elements.
 * This class assembles the weak form of:
 * \f$\mathbf{u} \cdot \nabla T - D \nabla^2T =S \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerDGCore : public TracerAssemblerBase<dim>
{
public:
  /**
   * @brief Default destructor.
   */
  virtual ~TracerAssemblerDGCore() = default;

  /**
   * @brief Assembles the matrix
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;


  /**
   * @brief Assembles the rhs
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;
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
  TracerAssemblerBDF(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Assembles the matrix
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */

  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;

  /**
   * @brief Assembles the rhs
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;

  // The simulation control is a necessary part of the transient terms.
  const std::shared_ptr<SimulationControl> simulation_control;
};

/**
 * @brief Assembles the symmetric interior penalty (SIPG) method (or
 * Nitsche's method) for internal faces. This assembler is only required
 * when solving the Tracer equation using a discontinuous Galerkin
 * discretization.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerSIPG : public TracerFaceAssemblerBase<dim>
{
public:
  TracerAssemblerSIPG()
  {}

  /**
   * @brief Interface for the call to matrix assembly
   * @param[in]  scratch_data Scratch data containing the Tracer
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */
  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedDGMethodsCopyData  &copy_data) override;


  /**
   * @brief Interface for the call to rhs
   * @param[in]  scratch_data Scratch data containing the Tracer
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedDGMethodsCopyData  &copy_data) override;
};

/**
 * @brief Assembles Nitsche's method for boundary faces.
 * This assembler is only required when solving the Tracer equation using a
 * discontinuous Galerkin discretization since the boundary conditions have
 * to be weakly imposed.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerBoundaryNitsche : public TracerFaceAssemblerBase<dim>
{
public:
  TracerAssemblerBoundaryNitsche(
    const BoundaryConditions::TracerBoundaryConditions<dim>
      &p_boundary_conditions_tracer)
    : boundary_conditions_tracer(p_boundary_conditions_tracer)
  {}

  /**
   * @brief Interface for the call to matrix assembly
   * @param[in]  scratch_data Scratch data containing the Tracer
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */
  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedDGMethodsCopyData  &copy_data) override;


  /**
   * @brief  Interface for the call to rhs
   * @param[in]  scratch_data Scratch data containing the Tracer
   * information. It is important to note that the scratch data has to have been
   * re-inited before calling for matrix assembly.
   * @param[in,out]  copy_data Destination where the local_rhs and local_matrix
   * should be copied
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedDGMethodsCopyData  &copy_data) override;

  BoundaryConditions::TracerBoundaryConditions<dim> boundary_conditions_tracer;
};

/**
 * @brief Class that assembles the reaction source term of the Tracer equation.
 * It assembles the term:
 * \f$ S = - R \= - k C \f$
 * Where C is the tracer concentration and k is a reaction rate prefactor that
 * changes depending on the underlying model. For power law of base rate
 * $\alpha$ and order $n$, the model is: \f$ S = - R \= - (\alpha C^(n-1)) C \f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class TracerAssemblerReaction : public TracerAssemblerBase<dim>
{
public:
  TracerAssemblerReaction(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief Assembles the matrix
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_matrix(const TracerScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;

  /**
   * @brief Assembles the rhs
   * @param[in] scratch_data (see base class)
   * @param[in,out] copy_data (see base class)
   */
  virtual void
  assemble_rhs(const TracerScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};



#endif
