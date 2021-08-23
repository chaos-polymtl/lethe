#include <dem/grid_motion.h>

#include <deal.II/grid/grid_tools.h>

#include <boost/math/special_functions.hpp>
#include <boost/range/adaptor/map.hpp>

using namespace dealii;

template <int dim>
GridMotion<dim>::GridMotion(
  const DEMSolverParameters<dim> &     dem_parameters,
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
           Parameters::Lagrangian::GridMotion<dim>::MotionType::free)
    {
      grid_motion = &GridMotion<dim>::move_grid_due_particles_forces;
      triangulation_inertia =
        dem_parameters.forces_torques.triangulation_inertia;
      boundary_rotational_velocity =
        dem_parameters.forces_torques.boundary_initial_rotational_velocity;
      boundary_translational_velocity =
        dem_parameters.forces_torques.boundary_initial_translational_velocity;
      GridMotion<dim>::pw_contact_force_object = pw_contact_force_object;
    }
}

template <>
void GridMotion<2>::move_grid_rotational(
  parallel::distributed::Triangulation<2> &triangulation)
{
  unsigned int count(0);
  for (unsigned int i = 0; i < 3; i++)
    {
      if (rotation_angle[i] != 0)
        {
          count++;
          if (count > 1)
            throw("In 2D, you cannot rotate around more than two axis");
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

template <int dim>
void
GridMotion<dim>::move_grid_due_particles_forces(
  parallel::distributed::Triangulation<dim> &triangulation)
{
  update_parameters_before_motion();
  calculate_motion_parameters();
  move_grid_translational(triangulation);
  // TODO : insert a move_grid_rotational(triangulation); The problem right now
  //  is that this function make the rotation around (0,0,0) point and not
  //  around the center of mass, which is moving
  update_parameters_after_motion();
}

template <int dim>
void
GridMotion<dim>::calculate_motion_parameters()
{
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));

  if (this_mpi_process == 0)
    {
      Tensor<1, dim> translational_velocity_one_step_time_further;
      Tensor<1, 3>   rotational_velocity_one_step_time_further;

      translational_velocity_one_step_time_further =
        (dt / triangulation_mass) * triangulation_forces +
        boundary_translational_velocity;

      for (unsigned int i = 0; i < (2 * dim - 3); i++)
        {
          rotational_velocity_one_step_time_further[i] =
            (dt / triangulation_inertia[i]) * triangulation_torques[i] +
            boundary_rotational_velocity[i];
        }

      // Shift vector should be the integral between two time step, here it's
      // done by taking the area (trapeze) formed by the velocity distance in a
      // time step.
      shift_vector = (translational_velocity_one_step_time_further +
                      boundary_translational_velocity) *
                     dt / 2;

      rotation_angle = (rotational_velocity_one_step_time_further +
                        boundary_rotational_velocity) *
                       dt / 2;

      boundary_translational_velocity =
        translational_velocity_one_step_time_further;
      boundary_rotational_velocity = rotational_velocity_one_step_time_further;
    }

  shift_vector   = Utilities::MPI::broadcast(MPI_COMM_WORLD, shift_vector);
  rotation_angle = Utilities::MPI::broadcast(MPI_COMM_WORLD, rotation_angle);
}

template <int dim>
void
GridMotion<dim>::update_parameters_before_motion()
{
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));

  if (this_mpi_process == 0)
    {
      // Update forces like torque and force
      std::map<unsigned int, Tensor<1, dim>> force_on_walls, torque_on_walls;
      force_on_walls  = pw_contact_force_object->get_force();
      torque_on_walls = pw_contact_force_object->get_torque();

      triangulation_forces  = 0;
      triangulation_torques = 0;
      for (auto it : force_on_walls)
        {
          triangulation_forces += it.second;
          triangulation_torques += torque_on_walls[it.first];
        }
    }
}

template <int dim>
void
GridMotion<dim>::update_parameters_after_motion()
{
  pw_contact_force_object->update_center_of_mass(shift_vector);
  // TODO : make an update_inertia(rotation_angle) function;
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
  for (auto &&pairs_in_contact_content :
       pw_pairs_in_contact | boost::adaptors::map_values)
    {
      for (auto &&contact_information :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          contact_information.normal_vector =
            updated_boundary_points_and_normal_vectors
              .at(contact_information.global_face_id)
              .first;
          contact_information.point_on_boundary =
            updated_boundary_points_and_normal_vectors
              .at(contact_information.global_face_id)
              .second;
        }
    }
}

template class GridMotion<2>;
template class GridMotion<3>;
