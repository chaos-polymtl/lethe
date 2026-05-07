// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/solid_constraints.h>

#include <deal.II/base/index_set.h>

#include <type_traits>
#include <unordered_set>

using namespace dealii;

namespace
{
  /**
   * @brief Flag DOFs connected to fluid cells.
   */
  template <int dim>
  void
  flag_dofs_connected_to_fluid(
    const FiniteElement<dim>                    &fe,
    const std::vector<types::global_dof_index>  &local_dof_indices,
    std::unordered_set<types::global_dof_index> &dofs_are_connected_to_fluid)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component = fe.system_to_component_index(i).first;
        if (component == dim)
          {
            dofs_are_connected_to_fluid.insert(local_dof_indices[i]);
          }
      }
  }

  /**
   * @brief Check if the cell is connected to a fluid cell.
   */
  bool
  check_cell_is_connected_to_fluid(
    const std::unordered_set<types::global_dof_index>
                                               &dofs_are_connected_to_fluid,
    const std::vector<types::global_dof_index> &local_dof_indices)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        auto search = dofs_are_connected_to_fluid.find(local_dof_indices[i]);
        if (search != dofs_are_connected_to_fluid.end())
          return true;
      }
    return false;
  }

  /**
   * @brief Constrain pressure DOFs if cells are not connected to fluid and
   * DOFs are locally owned.
   */
  template <int dim, typename DofsType>
  void
  constrain_pressure(
    const FiniteElement<dim>                   &fe,
    const bool                                 &non_zero_constraints,
    const DofsType                             &locally_owned_dofs,
    const std::vector<types::global_dof_index> &local_dof_indices,
    AffineConstraints<double>                  &constraints)
  {
    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      {
        const unsigned int component = fe.system_to_component_index(i).first;

        // Only pressure DOFs have an additional Dirichlet condition
        if (component == dim) // pressure DOFs
          {
            // We only apply the constraint on the locally owned pressure DOFs
            // since we have no way of verifying if the locally relevant DOFs
            // are connected to a fluid cell.
            bool dof_is_locally_owned = false;

            // For the GLS-family of solvers, we only have a single index set
            if constexpr (std::is_same_v<DofsType, IndexSet>)
              {
                dof_is_locally_owned =
                  locally_owned_dofs.is_element(local_dof_indices[i]);
              }

            // For the GD-family of solvers, we have two index sets. One for
            // velocities and one for pressure.
            if constexpr (std::is_same_v<DofsType, std::vector<IndexSet>>)
              {
                dof_is_locally_owned =
                  locally_owned_dofs[1].is_element(local_dof_indices[i]);
              }

            if (dof_is_locally_owned)
              {
                if (non_zero_constraints)
                  {
                    constraints.add_line(local_dof_indices[i]);
                    constraints.set_inhomogeneity(local_dof_indices[i], 0);
                  }
                else
                  {
                    constraints.add_line(local_dof_indices[i]);
                  }
              }
          }
      }
  }
} // namespace

template <int dim>
void
constrain_solid_cell_velocity_dofs(
  const FiniteElement<dim>                   &fe,
  const bool                                 &non_zero_constraints,
  const std::vector<types::global_dof_index> &local_dof_indices,
  AffineConstraints<double>                  &constraints)
{
  for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
    {
      const unsigned int component = fe.system_to_component_index(i).first;
      if (component < dim) // velocity DOFs
        {
          // We apply a constraint to all DOFs in the solid region, whether
          // they are locally owned or not.
          if (non_zero_constraints)
            {
              constraints.add_line(local_dof_indices[i]);
              constraints.set_inhomogeneity(local_dof_indices[i], 0);
            }
          else
            constraints.add_line(local_dof_indices[i]);
        }
    }
}

template <int dim, typename DofsType>
void
establish_solid_domain(const DoFHandler<dim>     &dof_handler,
                       const DofsType            &locally_owned_dofs,
                       const bool                 non_zero_constraints,
                       AffineConstraints<double> &constraints)
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // We will need to identify which pressure degrees of freedom are connected to
  // fluid region. For these, we won't establish a zero pressure constraint.
  std::unordered_set<types::global_dof_index> dofs_are_connected_to_fluid;

  // Loop through all cells to identify which cells are solid. This first step
  // is used to 1) constraint the velocity degree of freedom to be zero in the
  // solid region and 2) to identify which pressure degrees of freedom are
  // connected to fluid cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          // If the material_id is higher than 0, the region is a solid region.
          // Constrain the velocity DOFs to be zero.
          if (cell->material_id() > 0)
            {
              constrain_solid_cell_velocity_dofs(dof_handler.get_fe(),
                                                 non_zero_constraints,
                                                 local_dof_indices,
                                                 constraints);
            }
          else
            {
              // Cell is a fluid cell and as such all the pressure DOFs of that
              // cell are connected to the fluid. This will be used later on to
              // identify which pressure cells to constrain.
              flag_dofs_connected_to_fluid(dof_handler.get_fe(),
                                           local_dof_indices,
                                           dofs_are_connected_to_fluid);
            }
        }
    }

  // All pressure DOFs that are not connected to a fluid cell are constrained
  // to ensure that the system matrix has adequate conditioning.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          // If the material_id is > 0, the region is a solid region
          if (cell->material_id() > 0)
            {
              // First check if the cell is connected to a fluid cell by
              // checking if one of the DOF of the cell is connected to a fluid
              // cell.
              bool connected_to_fluid =
                check_cell_is_connected_to_fluid(dofs_are_connected_to_fluid,
                                                 local_dof_indices);

              // All the pressure DOFs with the cell are not connected to the
              // fluid. Consequently, we fix a constraint on these pressure
              // DOFs.
              if (!connected_to_fluid)
                {
                  constrain_pressure(dof_handler.get_fe(),
                                     non_zero_constraints,
                                     locally_owned_dofs,
                                     local_dof_indices,
                                     constraints);
                }
            }
        }
    }
}

// Explicit template instantiations
template void
constrain_solid_cell_velocity_dofs<2>(
  const FiniteElement<2>                     &fe,
  const bool                                 &non_zero_constraints,
  const std::vector<types::global_dof_index> &local_dof_indices,
  AffineConstraints<double>                  &constraints);
template void
constrain_solid_cell_velocity_dofs<3>(
  const FiniteElement<3>                     &fe,
  const bool                                 &non_zero_constraints,
  const std::vector<types::global_dof_index> &local_dof_indices,
  AffineConstraints<double>                  &constraints);

template void
establish_solid_domain<2, IndexSet>(const DoFHandler<2>       &dof_handler,
                                    const IndexSet            &locally_owned_dofs,
                                    const bool                 non_zero_constraints,
                                    AffineConstraints<double> &constraints);
template void
establish_solid_domain<3, IndexSet>(const DoFHandler<3>       &dof_handler,
                                    const IndexSet            &locally_owned_dofs,
                                    const bool                 non_zero_constraints,
                                    AffineConstraints<double> &constraints);
template void
establish_solid_domain<2, std::vector<IndexSet>>(
  const DoFHandler<2>          &dof_handler,
  const std::vector<IndexSet>  &locally_owned_dofs,
  const bool                    non_zero_constraints,
  AffineConstraints<double>    &constraints);
template void
establish_solid_domain<3, std::vector<IndexSet>>(
  const DoFHandler<3>          &dof_handler,
  const std::vector<IndexSet>  &locally_owned_dofs,
  const bool                    non_zero_constraints,
  AffineConstraints<double>    &constraints);
