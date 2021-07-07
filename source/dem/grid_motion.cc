#include <dem/grid_motion.h>

#include <deal.II/grid/grid_tools.h>

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
  else if (dem_parameters.grid_motion.motion_type ==
           Parameters::Lagrangian::GridMotion<
             dim>::MotionType::translational_rotational)
    {
      grid_motion = &GridMotion<dim>::move_grid_translational_rotational;
      rotation_angle =
        dem_parameters.grid_motion.grid_rotational_speed * dem_time_step;
      rotation_axis = dem_parameters.grid_motion.grid_rotational_axis;
      shift_vector =
        dem_parameters.grid_motion.grid_translational_velocity * dem_time_step;
    }
  else
    {
      throw std::runtime_error("Specified grid motion is not valid");
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

template <>
void GridMotion<2>::move_grid_translational_rotational(
  parallel::distributed::Triangulation<2> &triangulation)
{
  GridTools::shift(shift_vector, triangulation);
  GridTools::rotate(rotation_angle, triangulation);
}

template <>
void GridMotion<3>::move_grid_translational_rotational(
  parallel::distributed::Triangulation<3> &triangulation)
{
  GridTools::shift(shift_vector, triangulation);
  GridTools::rotate(rotation_angle, rotation_axis, triangulation);
}

template class GridMotion<2>;
template class GridMotion<3>;
