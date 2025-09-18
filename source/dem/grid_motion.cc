// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/grid_motion.h>

#include <deal.II/grid/grid_tools.h>

using namespace dealii;

template <int dim, int spacedim>
GridMotion<dim, spacedim>::GridMotion(
  const Parameters::Lagrangian::GridMotion<spacedim> &grid_motion_parameters,
  const double                                        dem_time_step)
{
  auto *action_manager = DEMActionManager::get_action_manager();

  switch (grid_motion_parameters.motion_type)
    {
      case Parameters::Lagrangian::GridMotion<spacedim>::MotionType::none:
        {
          grid_motion = &GridMotion<dim, spacedim>::no_motion;
          break;
        }
      case Parameters::Lagrangian::GridMotion<spacedim>::MotionType::rotational:
        {
          // Communicate to the action manager that there is a grid motion
          action_manager->set_grid_motion_enabled();

          grid_motion = &GridMotion<dim, spacedim>::move_grid_rotational;
          rotation_angle =
            grid_motion_parameters.grid_rotational_speed * dem_time_step;
          rotation_axis = grid_motion_parameters.grid_rotational_axis;
          break;
        }
      case Parameters::Lagrangian::GridMotion<
        spacedim>::MotionType::translational:
        {
          // Flag and communicate to the action manager that there is a
          // grid motion
          action_manager->set_grid_motion_enabled();

          grid_motion = &GridMotion<dim, spacedim>::move_grid_translational;
          shift_vector =
            grid_motion_parameters.grid_translational_velocity * dem_time_step;
          break;
        }
      default:
        AssertThrow(false, ExcNotImplemented());
    }
}

template <>
void
GridMotion<1, 2>::move_grid_rotational(Triangulation<1, 2> &)
{
  throw ExcImpossibleInDim(1);
  // TODO We need to add this function to GridTools for dim=1
  // GridTools::rotate(rotation_angle, triangulation);
}

template <>
void
GridMotion<2, 2>::move_grid_rotational(Triangulation<2, 2> &triangulation)
{
  GridTools::rotate(rotation_angle, triangulation);
}

template <>
void
GridMotion<2, 3>::move_grid_rotational(Triangulation<2, 3> &triangulation)
{
  Tensor<1, 3> axis({0, 0, 0});
  axis[rotation_axis] = 1;
  GridTools::rotate(axis, rotation_angle, triangulation);
}

template <>
void
GridMotion<3, 3>::move_grid_rotational(Triangulation<3, 3> &triangulation)
{
  Tensor<1, 3> axis({0, 0, 0});
  axis[rotation_axis] = 1;
  GridTools::rotate(axis, rotation_angle, triangulation);
}

template <int dim, int spacedim>
void
GridMotion<dim, spacedim>::move_grid_translational(
  Triangulation<dim, spacedim> &triangulation)
{
  GridTools::shift(shift_vector, triangulation);
}

template <int dim, int spacedim>
void
GridMotion<dim, spacedim>::
  update_boundary_points_and_normal_vectors_in_contact_list(
    typename DEM::dem_data_structures<spacedim>::particle_wall_in_contact
      &particle_wall_pairs_in_contact,
    const typename DEM::dem_data_structures<
      spacedim>::boundary_points_and_normal_vectors
      &updated_boundary_points_and_normal_vectors)
{
  // If there is no grid motion, exit the function
  if (!DEMActionManager::get_action_manager()->check_grid_motion_enabled())
    return;

  for (auto &[particle_id, pairs_in_contact_content] :
       particle_wall_pairs_in_contact)
    {
      // Prevent compiler warning
      (void)particle_id;
      for (auto pairs_in_contact_iterator = pairs_in_contact_content.begin();
           pairs_in_contact_iterator != pairs_in_contact_content.end();)
        {
          auto  global_face_id      = pairs_in_contact_iterator->first;
          auto &contact_information = pairs_in_contact_iterator->second;

          // Since we used the negative keys for diamond-shaped cells, we check
          // and remove these elements from the particle_wall_pairs_in_contact
          if (global_face_id >= 0)
            {
              contact_information.normal_vector =
                updated_boundary_points_and_normal_vectors.at(global_face_id)
                  .first;
              contact_information.point_on_boundary =
                updated_boundary_points_and_normal_vectors.at(global_face_id)
                  .second;
              ++pairs_in_contact_iterator;
            }
          else
            {
              pairs_in_contact_content.erase(pairs_in_contact_iterator++);
            }
        }
    }
}

template class GridMotion<1, 2>;
template class GridMotion<2, 2>;
template class GridMotion<2, 3>;
template class GridMotion<3, 3>;
