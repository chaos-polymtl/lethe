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
 * @param[in] local_dof_indices Vector of a cell's local DOF indices.
 *
 * @param[in] non_zero_constraints If this parameter is true, it indicates
 * that the constraints used are non-zero constraints. If
 * this is set to false, the constraints are zero constraints.
 *
 * @param[out] constraints Homogeneous constraints holding object.
 */
template <int dim>
void
constrain_solid_cell_velocity_dofs(
  const dealii::FiniteElement<dim>                   &fe,
  const std::vector<dealii::types::global_dof_index> &local_dof_indices,
  const bool                                         &non_zero_constraints,
  dealii::AffineConstraints<double>                  &constraints);

/**
 * @brief Turn regions of the mesh where the @p material_id>0 into a solid
 * block by injecting velocity and pressure DOFs into the constraints.
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
 *
 * @param[in] locally_owned_dofs The locally owned DoFs for the dof_handler.
 * Used to filter pressure DoFs that should receive a
 * Dirichlet constraint when their cell is fully
 * disconnected from the fluid.
 *
 * @param[in] non_zero_constraints If this parameter is true, it indicates
 * that the constraints used for the solid domainare the non-zero constraints.
 * If this is set to false, the constraints being used are zero constraints.
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

/**
 * @brief Multigrid-level variant of establish_solid_domain for the local
 * smoothing multigrid (LSMG) preconditioner.
 *
 * Unlike the active-cell version above, this overload walks the cells of a
 * specific multigrid level using @p mg_cell_iterators_on_level and uses
 * @p get_mg_dof_indices, so the DOF indices it adds to the constraint object
 * match the multigrid-level IndexSet with which @p constraints was reinit'd
 * (i.e. @p locally_owned_mg_dofs(level) and the corresponding locally
 * relevant level DoFs). Using the active-cell version on an LSMG level
 * constraint object triggers the @p local_lines.is_element(line_n) assertion
 * inside @p AffineConstraints::add_line because active and MG DoF numberings
 * are unrelated.
 *
 * The @p material_id is inherited from active cells by their parents on
 * coarser levels, so the solid/fluid classification is well defined on MG
 * levels.
 *
 * @param[in] dof_handler DoFHandler distributed with @p distribute_mg_dofs.
 *
 * @param[in] level Multigrid level on which constraints are built.
 *
 * @param[in] locally_owned_mg_dofs The locally owned MG DoFs on @p level
 * (i.e. @p dof_handler.locally_owned_mg_dofs(level)). Used to filter pressure
 * DoFs that should receive a Dirichlet constraint when their cell is fully
 * disconnected from the fluid.
 *
 * @param[in] non_zero_constraints If true, set inhomogeneous (zero-valued)
 * constraints; otherwise add homogeneous constraints.
 *
 * @param[in,out] constraints Level constraint object to populate.
 */
template <int dim>
void
establish_solid_domain_lsmg(const dealii::DoFHandler<dim> &dof_handler,
                            const unsigned int             level,
                            const dealii::IndexSet &locally_owned_mg_dofs,
                            const bool              non_zero_constraints,
                            dealii::AffineConstraints<double> &constraints);

#endif
