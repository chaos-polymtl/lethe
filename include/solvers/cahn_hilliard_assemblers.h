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

#ifndef lethe_cahn_hilliard_assemblers_h
#define lethe_cahn_hilliard_assemblers_h

#include <core/simulation_control.h>

#include <solvers/cahn_hilliard_scratch_data.h>
#include <solvers/copy_data.h>
#include <solvers/multiphysics_interface.h>

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
    const std::shared_ptr<SimulationControl> simulation_control)
    : simulation_control(simulation_control)
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
  const std::shared_ptr<SimulationControl> simulation_control;
};


/**
 * @brief Class that assembles the core of the Cahn-Hilliard equation :
 * \f$ \frac{d \phi}{dt} +  (u \cdot \nabla) \phi - \nabla \cdot (M(\phi)\nabla
 * \eta) = 0 \\ \eta -  \frac{\lambda}{\epsilon^2}(\phi^3 - \phi) + \lambda
 * \nabla^2 \phi + \xi \nabla^2 \eta  = 0 \f$
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
    const std::shared_ptr<SimulationControl> simulation_control,
    const Parameters::CahnHilliard           cahn_hilliard_parameters,
    const double                             epsilon)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , cahn_hilliard_parameters(cahn_hilliard_parameters)
    , epsilon(epsilon)
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

  const Parameters::CahnHilliard cahn_hilliard_parameters;
  // Epsilon is a coefficient which depends on the mesh size. The thickness of
  // the interface between the two phases is proportionnal to espilon
  const double epsilon;
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
    const std::shared_ptr<SimulationControl> simulation_control,
    const Parameters::CahnHilliard           cahn_hilliard_parameters,
    const double                             epsilon,
    const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
      &p_boundary_conditions_cahn_hilliard)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , cahn_hilliard_parameters(cahn_hilliard_parameters)
    , epsilon(epsilon)
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

  const Parameters::CahnHilliard cahn_hilliard_parameters;
  // Epsilon is a coefficient which depends on the mesh size. The thickness of
  // the interface between the two phases is proportionnal to espilon
  const double epsilon;
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
    const std::shared_ptr<SimulationControl> simulation_control,
    const Parameters::CahnHilliard           cahn_hilliard_parameters,
    const double                             epsilon,
    const BoundaryConditions::CahnHilliardBoundaryConditions<dim>
      &p_boundary_conditions_cahn_hilliard)
    : CahnHilliardAssemblerBase<dim>(simulation_control)
    , cahn_hilliard_parameters(cahn_hilliard_parameters)
    , epsilon(epsilon)
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


  const Parameters::CahnHilliard cahn_hilliard_parameters;
  // Epsilon is a coefficient which depends on the mesh size. The thickness of
  // the interface between the two phases is proportionnal to espilon
  const double epsilon;
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
    const std::shared_ptr<SimulationControl> simulation_control)
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
