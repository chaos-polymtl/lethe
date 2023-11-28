#include <dem/find_maximum_particle_size.h>
#include <dem/input_parameter_inspection.h>

using namespace dealii;

template <int dim>
void
input_parameter_inspection(
  const DEMSolverParameters<dim>              &dem_parameters,
  const ConditionalOStream                    &pcout,
  const std::vector<shared_ptr<Distribution>> &distribution_object_container)
{
  // Getting the input parameters as local variable
  auto   parameters          = dem_parameters;
  auto   physical_properties = dem_parameters.lagrangian_physical_properties;
  double rayleigh_time_step  = 1. / DBL_MIN;
  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      double shear_modulus =
        physical_properties.youngs_modulus_particle[i] /
        (2.0 * (1.0 + physical_properties.poisson_ratio_particle[i]));

      double min_diameter =
        distribution_object_container.at(i)->find_min_diameter();

      rayleigh_time_step = std::min(
        M_PI_2 * min_diameter *
          sqrt(physical_properties.density_particle[i] / shear_modulus) /
          (0.1631 * physical_properties.poisson_ratio_particle[i] + 0.8766),
        rayleigh_time_step);
    }

  const double time_step_rayleigh_ratio =
    parameters.simulation_control.dt / rayleigh_time_step;
  pcout << "DEM time-step is " << time_step_rayleigh_ratio * 100
        << "% of Rayleigh time step" << std::endl;

  if (time_step_rayleigh_ratio > 0.15)
    {
      pcout << "Warning: It is recommended to decrease the time-step"
            << std::endl;
    }
  else if (time_step_rayleigh_ratio < 0.01)
    {
      pcout << "Warning: It is recommended to increase the time-step"
            << std::endl;
    }

  // Checking particle size range
  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    {
      if (distribution_object_container.at(i)->find_min_diameter() < 0.)
        {
          pcout
            << "Warning: Requested particle size distribution for type: " << i
            << " is not well-defined. Using requested distribution may lead to "
               "changes in sampled particle sizes to avoid negative particle "
               "diameters. You can consider decreasing the standard deviation of "
               "size distribution or increasing the particle diameter instead."
            << std::endl;
        }
    }

  // Insertion parameters check
  const double insertion_distance_per_particle =
    0.5 * (parameters.insertion_info.distance_threshold - 1);

  if (parameters.insertion_info.insertion_method ==
        Parameters::Lagrangian::InsertionInfo::InsertionMethod::volume &&
      parameters.insertion_info.random_number_range >=
        insertion_distance_per_particle)
    pcout
      << "Warning: Particles may have collision at the insertion step. This can lead"
         " to high initial velocities (due to initial overlap) or errors when using "
         "less stable integration schemes. It is recommended to decrease the random "
         "number range or to increase the insertion distance threshold."
      << std::endl;

  // Grid motion check
  if (parameters.grid_motion.motion_type ==
      Parameters::Lagrangian::GridMotion<dim>::MotionType::rotational)
    {
      if (dim == 2)
        {
          if (parameters.grid_motion.grid_rotational_axis != 0 &&
              parameters.grid_motion.grid_rotational_axis != 1)
            throw std::runtime_error(
              "Specified grid rotational axis is not valid, use 0 for rotation around x"
              " axis or 1 for rotation around y axis.");
        }
      else if (dim == 3)
        {
          if (parameters.grid_motion.grid_rotational_axis != 0 &&
              parameters.grid_motion.grid_rotational_axis != 1 &&
              parameters.grid_motion.grid_rotational_axis != 2)
            throw std::runtime_error(
              "Specified grid rotational axis is not valid, use 0 for rotation around x"
              " axis or 1 for rotation around y axis or 2 for rotation around z axis.");
        }
    }
}

template void
input_parameter_inspection(
  const DEMSolverParameters<2>                &dem_parameters,
  const ConditionalOStream                    &pcout,
  const std::vector<shared_ptr<Distribution>> &distribution_object_container);

template void
input_parameter_inspection(
  const DEMSolverParameters<3>                &dem_parameters,
  const ConditionalOStream                    &pcout,
  const std::vector<shared_ptr<Distribution>> &distribution_object_container);
