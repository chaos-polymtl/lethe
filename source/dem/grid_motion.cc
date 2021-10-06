#include <dem/grid_motion.h>

#include <deal.II/grid/grid_tools.h>

#include <boost/math/special_functions.hpp>
#include <boost/range/adaptor/map.hpp>

using namespace dealii;

template <int dim>
GridMotion<dim>::GridMotion(const DEMSolverParameters<dim> &dem_parameters,
                            const double &                  dem_time_step)
{
  // Setting grid motion type
  if (dem_parameters.grid_motion.motion_type ==
      Parameters::Lagrangian::GridMotion<dim>::MotionType::rotational)
    {
      grid_motion = &GridMotion<dim>::move_grid_rotational;
      rotation_angle =
        dem_parameters.grid_motion.grid_rotational_speed * dem_time_step;
      rotation_axis = dem_parameters.grid_motion.grid_rotational_axis;
    }
  else if (dem_parameters.grid_motion.motion_type ==
           Parameters::Lagrangian::GridMotion<dim>::MotionType::translational)
    {
      grid_motion = &GridMotion<dim>::move_grid_translational;
      shift_vector =
        dem_parameters.grid_motion.grid_translational_velocity * dem_time_step;
    }
}

template <>
void GridMotion<2>::move_grid_rotational(
  parallel::distributed::Triangulation<2> &triangulation)
{
  GridTools::rotate(rotation_angle, triangulation);
}

template <>
void GridMotion<3>::move_grid_rotational(
  parallel::distributed::Triangulation<3> &triangulation)
{
  GridTools::rotate(rotation_angle, rotation_axis, triangulation);
}

template <int dim>
void
GridMotion<dim>::move_grid_translational(
  parallel::distributed::Triangulation<dim> &triangulation)
{
  GridTools::shift(shift_vector, triangulation);
}

template <int dim>
void
GridMotion<dim>::update_boundary_points_and_normal_vectors_in_contact_list(
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    &pw_pairs_in_contact,
  const std::map<unsigned int, std::pair<Tensor<1, dim>, Point<dim>>>
    &updated_boundary_points_and_normal_vectors)
{
  for (auto &[particle_id, pairs_in_contact_content] : pw_pairs_in_contact)
    {
      for (auto pairs_in_contact_iterator = pairs_in_contact_content.begin();
           pairs_in_contact_iterator != pairs_in_contact_content.end();)
        {
          auto &contact_information = pairs_in_contact_iterator->second;

          // Since we used the negative keys for diamond-shaped cells, we check
          // and remove these elements from the pw_pairs_in_contact
          if (contact_information.global_face_id >= 0)
            {
              contact_information.normal_vector =
                updated_boundary_points_and_normal_vectors
                  .at(contact_information.global_face_id)
                  .first;
              contact_information.point_on_boundary =
                updated_boundary_points_and_normal_vectors
                  .at(contact_information.global_face_id)
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

template class GridMotion<2>;
template class GridMotion<3>;
