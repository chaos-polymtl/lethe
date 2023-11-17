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

#include <solvers/cahn_hilliard_scratch_data.h>
#include <solvers/copy_data.h>
#include <solvers/multiphysics_interface.h>


#ifndef lethe_cahn_hilliard_assemblers_h
#  define lethe_cahn_hilliard_assemblers_h

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Cahn-Hilliard equation
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CahnHilliardAssemblerBase
{
public:
  CahnHilliardAssemblerBase(
    std::shared_ptr<SimulationControl> p_simulation_control)
    : simulation_control(p_simulation_control)
  {}

  /**
   * @brief assemble_matrix Interface for the call to matrix assembly
   * @param scratch_data Scratch data containing the information for Cahn-Hilliard
   * equations.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for matrix assembly.
   * @param copy_data Destination where the local_rhs and local_matrix are copied to.
   */

  virtual void
  assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) = 0;


  /**
   * @brief assemble_rhs Interface for the call to rhs assembly
   * @param scratch_data Scratch data containing the information for the Cahn-Hilliard
   * equations.
   * It is important to note that the scratch data has to have been re-inited
   * before calling for rhs assembly.
   * @param copy_data Destination where the local_rhs and local_matrix are copied to.
   */

  virtual void
  assemble_rhs(CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) = 0;

protected:
  std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the core of the Cahn-Hilliard equation
 * with the following weak form:
 * dPhi/dt +  u * gradPhi =  div(M(Phi)*grad eta)
 * eta - f(Phi) + epsilon^2 * div(grad Phi) = 0
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */


template <int dim>
class CahnHilliardAssemblerCore : public CahnHilliardAssemblerBase<dim>
{
public:
  CahnHilliardAssemblerCore(
    std::shared_ptr<SimulationControl> simulation_control,
    Parameters::CahnHilliard           cahn_hilliard_parameters)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , cahn_hilliard_parameters(cahn_hilliard_parameters)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;

  Parameters::CahnHilliard cahn_hilliard_parameters;
};


/**
 * @brief Class that assembles the boundary condition on the angle of contact in the
 * Cahn-Hilliard equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CahnHilliardAssemblerAngleOfContact
  : public CahnHilliardAssemblerBase<dim>
{
public:
  CahnHilliardAssemblerAngleOfContact(
    std::shared_ptr<SimulationControl> simulation_control,
    const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
      &p_boundary_conditions_cahn_hilliard)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , boundary_conditions_cahn_hilliard(p_boundary_conditions_cahn_hilliard)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;


  Parameters::MaterialInteractions material_interaction_parameters;
  const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
    &boundary_conditions_cahn_hilliard;
};

/**
 * @brief Class that assembles the boundary condition that allows a free angle of
 * contact, thus adding a degree of freedom to the problem
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CahnHilliardAssemblerFreeAngle : public CahnHilliardAssemblerBase<dim>
{
public:
  CahnHilliardAssemblerFreeAngle(
    std::shared_ptr<SimulationControl> simulation_control,
    const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
      &p_boundary_conditions_cahn_hilliard)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , boundary_conditions_cahn_hilliard(p_boundary_conditions_cahn_hilliard)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;


  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;


  const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
    &boundary_conditions_cahn_hilliard;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Cahn-Hilliard equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{T}^{t+\Delta t}-\mathbf{T}^{t}}{\Delta t}\f$
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class CahnHilliardAssemblerBDF : public CahnHilliardAssemblerBase<dim>
{
public:
  CahnHilliardAssemblerBDF(
    std::shared_ptr<SimulationControl> simulation_control)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData    &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData    &copy_data) override;
};


#endif
