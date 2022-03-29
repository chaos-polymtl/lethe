#include <dem/grid_motion.h>

#include <deal.II/grid/grid_tools.h>

#include <boost/math/special_functions.hpp>
#include <boost/range/adaptor/map.hpp>

#include <math.h>

using namespace dealii;

template <int dim>
GridMotion<dim>::GridMotion(
  DEMSolverParameters<dim> &           dem_parameters,
  const double &                       dem_time_step,
  std::shared_ptr<PWContactForce<dim>> pw_contact_force_object)
  : triangulation_mass(dem_parameters.forces_torques.triangulation_mass)
  , dt(dem_parameters.simulation_control.dt)
{
  // Setting grid motion type
  if (dem_parameters.grid_motion.motion_type ==
      Parameters::Lagrangian::GridMotion<dim>::MotionType::rotational)
    {
      grid_motion = &GridMotion<dim>::move_grid_rotational;
      rotation_angle =
        dem_parameters.grid_motion.grid_rotational_speed * dem_time_step;
    }
  else if (dem_parameters.grid_motion.motion_type ==
           Parameters::Lagrangian::GridMotion<dim>::MotionType::translational)
    {
      grid_motion = &GridMotion<dim>::move_grid_translational;
      shift_vector =
        dem_parameters.grid_motion.grid_translational_velocity * dem_time_step;
    }
  else if (dem_parameters.grid_motion.motion_type ==
           Parameters::Lagrangian::GridMotion<dim>::MotionType::cylinder_motion)
    {
      grid_motion = &GridMotion<dim>::cylinder_motion;
      triangulation_inertia =
        dem_parameters.forces_torques.triangulation_inertia;
      GridMotion<dim>::pw_contact_force_object = pw_contact_force_object;
      rotation_axis = dem_parameters.grid_motion.cylinder_rotation_axis;

      for (unsigned int d = 0; d < dim; ++d)
        {
          if (d == rotation_axis)
            cylinder_rotation_unit_vector[d] = 1;
          else
            cylinder_rotation_unit_vector[d] = 0;
        }

      // If the grid motion type is cylinder_motion, then calculate_force_torque
      // should be enabled
      dem_parameters.forces_torques.calculate_force_torque = true;

      for (unsigned int d = 0; d < dim; ++d)
        gravity[d] = dem_parameters.lagrangian_physical_properties.g[d];

      inclined_plane_angle =
        dem_parameters.grid_motion.inclined_plane_angle * M_PI / 180;

      cylinder_radius = dem_parameters.grid_motion.cylinder_radius;
    }
}

template <>
void GridMotion<2>::move_grid_rotational(
  parallel::distributed::Triangulation<2> &triangulation)
{
  unsigned int count(0);
  for (unsigned int i = 0; i < 3; i++)
    {
      if (rotation_angle[i] != 0.0)
        {
          count++;
          if (count > 1)
            throw std::runtime_error(
              "In 2D, you cannot rotate around more than two axis");
          GridTools::rotate(rotation_angle[i], triangulation);
        }
    }
}

template <>
void GridMotion<3>::move_grid_rotational(
  parallel::distributed::Triangulation<3> &triangulation)
{
  for (unsigned int i = 0; i < 3; i++)
    {
      GridTools::rotate(rotation_angle[i], i, triangulation);
    }
}

template <int dim>
void
GridMotion<dim>::move_grid_translational(
  parallel::distributed::Triangulation<dim> &triangulation)
{
  GridTools::shift(shift_vector, triangulation);
}

// Note that the cylinder_motion is only available for rotation of
// three-dimensional cylindrical geometries at the moment.
template <>
void GridMotion<3>::cylinder_motion(
  parallel::distributed::Triangulation<3> &triangulation)
{
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));

  // the colorization of the cylinder hull is 0, so we only use the forces and
  // torques acting on the hull
  const unsigned int cylinder_hull_colorization = 0;

  // Calculation of rotation angle is only performed on one process in
  // parallel simulations
  if (this_mpi_process == 0)
    {
      // Update torque and force on walls
      std::map<unsigned int, Tensor<1, 3>> force_on_walls, torque_on_walls;

      // Call get_force() get_torque() functions in PWContactForce class
      force_on_walls  = pw_contact_force_object->get_force();
      torque_on_walls = pw_contact_force_object->get_torque();

      triangulation_torques = 0;

      // The angular acceleration acting on the rolling cylinder is equal to
      // $$\alpha = (\mathbf{T_{particles}} - (\mathbf{F_{particles}} + m *
      // \mathbf{g}) * \sin(\theta) * R) / \mathbf{I_{eff}}$$ Since the gravity
      // is already multiplied to $$\sin(\thetha)$$ in the constructor of the
      // dem class, we only multiply the
      // $$\mathbf{F_{particles}}$$ to $$\sin(\thetha)$$. Then we multiply the
      // summation of the gravitational and $$\mathbf{F_{particles}}$$ to
      // cylinder_rotation_unit_vector.

      triangulation_torques = torque_on_walls[cylinder_hull_colorization] -
                              (force_on_walls[cylinder_hull_colorization][1] *
                                 sin(inclined_plane_angle) +
                               triangulation_mass * gravity.norm()) *
                                cylinder_rotation_unit_vector * cylinder_radius;


      Tensor<1, 3> rotational_velocity_one_time_step_further;


      for (unsigned int i = 0; i < 3; i++)
        {
          if (triangulation_inertia[i] != 0)
            rotational_velocity_one_time_step_further[i] =
              (dt / triangulation_inertia[i]) * triangulation_torques[i] +
              boundary_rotational_velocity[i];
          else
            throw std::runtime_error(
              "Triangulation inertia cannot be equal to zero");
        }

      rotation_angle = (rotational_velocity_one_time_step_further +
                        boundary_rotational_velocity) *
                       dt / 2;

      boundary_rotational_velocity = rotational_velocity_one_time_step_further;
    }

  // The following line enables locking the rotation of the cylinder around
  // desired axis (cylinder_rotation_unit_vector which is equivalent to
  // "cylinder rotation axis" parameter in the parameter handler.
  double rotation_angle_around_axis =
    rotation_angle * cylinder_rotation_unit_vector;

  // Broadcast the rotation angle to all the processes
  rotation_angle_around_axis =
    Utilities::MPI::broadcast(MPI_COMM_WORLD, rotation_angle_around_axis);

  // Apply the rotation using rotate function
  GridTools::rotate(rotation_angle_around_axis, rotation_axis, triangulation);
}

template <>
void GridMotion<2>::cylinder_motion(
  parallel::distributed::Triangulation<2> & /*triangulation*/)
{}

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
      // Prevent compiler warning
      (void)particle_id;
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
