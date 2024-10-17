// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cahn_hilliard_assemblers_h
#define lethe_cahn_hilliard_assemblers_h

#include <core/simulation_control.h>

#include <solvers/cahn_hilliard_scratch_data.h>
#include <solvers/copy_data.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/physics_assemblers.h>

/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Cahn-Hilliard equation.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup assemblers
 */
template <int dim>
using CahnHilliardAssemblerBase =
  PhysicsAssemblerBase<CahnHilliardScratchData<dim>, StabilizedMethodsCopyData>;

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
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::CahnHilliard            cahn_hilliard_parameters,
    const double                              epsilon)
    : simulation_control(simulation_control)
    , cahn_hilliard_parameters(cahn_hilliard_parameters)
    , epsilon(epsilon)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CahnHilliard           cahn_hilliard_parameters;
  // Epsilon is a coefficient that depends on the mesh size. The thickness of
  // the interface between the two phases is proportional to espilon
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
    : simulation_control(simulation_control)
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
  assemble_matrix(const CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CahnHilliard           cahn_hilliard_parameters;
  // Epsilon is a coefficient that depends on the mesh size. The thickness of
  // the interface between the two phases is proportional to espilon
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
    : simulation_control(simulation_control)
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
  assemble_matrix(const CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
  const Parameters::CahnHilliard           cahn_hilliard_parameters;
  // Epsilon is a coefficient that depends on the mesh size. The thickness of
  // the interface between the two phases is proportional to espilon
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
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const CahnHilliardScratchData<dim> &scratch_data,
                  StabilizedMethodsCopyData          &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const CahnHilliardScratchData<dim> &scratch_data,
               StabilizedMethodsCopyData          &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};


#endif
