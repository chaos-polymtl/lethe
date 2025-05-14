// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef rans_assemblers_h
#define rans_assemblers_h

#include <core/boundary_conditions.h>
#include <core/simulation_control.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_scratch_data.h>
#include <solvers/physics_assemblers.h>

/*
 * Exceptions used to capture incoherent setup of assemblers
 */



/**
 * @brief A pure virtual class that serves as an interface for all
 * of the assemblers for the Navier-Stokes Equations
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
using NavierStokesAssemblerBase =
  PhysicsAssemblerBase<NavierStokesScratchData<dim>,
                       StabilizedMethodsTensorCopyData<dim>>;


/**
 * @brief Class that assembles the core of the Navier-Stokes equation.
 * According to the following weak form:
 * \f$\mathbf{u} \cdot \nabla \mathbf{u} - \nabla p - \nu \nabla^2 \mathbf{u}
 * =0 \f$ with an SUPG and PSPG stabilization
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class PSPGSUPGRANSTurbulenceAssemblerCore
  : public NavierStokesAssemblerBase<dim>
{
public:
  PSPGSUPGRANSTurbulenceAssemblerCore(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};

/**
 * @brief Class that assembles the transient time arising from BDF time
 * integration for the Navier Stokes equations. For example, if a BDF1 scheme is
 * chosen, the following is assembled
 * \f$\frac{\mathbf{u}^{t+\Delta t}-\mathbf{u}^{t}{\Delta t}
 *
 * @tparam dim An integer that denotes the number of spatial dimensions
 *
 * @ingroup assemblers
 */
template <int dim>
class PSPGSUPGRANSTurbulenceAssemblerBDF : public NavierStokesAssemblerBase<dim>
{
public:
  PSPGSUPGRANSTurbulenceAssemblerBDF(
    const std::shared_ptr<SimulationControl> &simulation_control)
    : simulation_control(simulation_control)
  {}

  /**
   * @brief assemble_matrix Assembles the matrix
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */

  virtual void
  assemble_matrix(const NavierStokesScratchData<dim>   &scratch_data,
                  StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief assemble_rhs Assembles the rhs
   * @param scratch_data (see base class)
   * @param copy_data (see base class)
   */
  virtual void
  assemble_rhs(const NavierStokesScratchData<dim>   &scratch_data,
               StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  const std::shared_ptr<SimulationControl> simulation_control;
};

#endif
