// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_solid_constraints_h
#define lethe_solid_constraints_h

#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/affine_constraints.h>

#include <vector>

/**
 * @brief Constrain velocity DOFs of a solid cell.
 *
 * @param [in] fe Finite element used to extract the different components.
 * Inherently, we assume that this fe has dim+1 component and corresponds
 * to a Navier-Stokes problem.
 *
 * @param[in] non_zero_constraints If this parameter is true, it indicates
 * that non-zero constraints are applied in the solid domain. If
 * this is set to false, homogeneous constraints are applied in the solid
 * domain.
 *
 * @param[in] local_dof_indices Vector of a cell's local DOF indices.
 *
 * @param[out] constraints Homogeneous constraints holding object.
 */
template <int dim>
void
constrain_solid_cell_velocity_dofs(
  const dealii::FiniteElement<dim>                   &fe,
  const bool                                         &non_zero_constraints,
  const std::vector<dealii::types::global_dof_index> &local_dof_indices,
  dealii::AffineConstraints<double>                  &constraints);

/**
 * @brief Turn regions of the mesh where the @p material_id>0 into a solid
 * block by injecting velocity and pressure DOFs into the zero constraints.
 *
 * It is achieved by imposing \f$\mathbf{u}=0\f$ within the cells which have a
 * @p material_id>0. In addition, solid cells which are not connected to the
 * fluid by any means also get a pressure Dirichlet boundary condition which
 * fixes the pressure to 0. It ensures that the linear system is well-posed.
 * Right now, this routine only supports the usage of 1 solid domain, but
 * eventually it could be extended to more than one. By default, the fluid
 * domain is assumed to have a @p material_id=0 and the rest of the domains
 * have a @p material_id>0.
 *
 * @param[in] dof_handler DoFHandler for which the solid region is
 * established. In the majority of cases this is the DoFHandler that is a
 * member of this class, but if a multigrid preconditioner, constraints must
 * be built at every level. In this case, this DoFHandler can be the
 * DoFHandler at any of the multigrid level.
 *
 * @param[in] non_zero_constraints If this parameter is true, it indicates
 * that non-zero constraints are being constrained for the solid domain. If
 * this is set to false, homogeneous constraints are constrained in the solid
 * domain.
 *
 * @param[in] constraints Set of constraints on which the solid domain
 * constraint is applied.
 */
template <int dim, typename DofsType>
void
establish_solid_domain(const dealii::DoFHandler<dim>     &dof_handler,
                       const DofsType                    &locally_owned_dofs,
                       const bool                         non_zero_constraints,
                       dealii::AffineConstraints<double> &constraints);

#endif
