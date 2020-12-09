#include <dem/find_maximum_particle_size.h>
#include <dem/input_parameter_inspection.h>

using namespace dealii;

template <int dim>
void
input_parameter_inspection(const DEMSolverParameters<dim> &dem_parameters,
                           const ConditionalOStream &      pcout,
                           const double &standard_deviation_multiplier)
{
  // Getting the input parameters as local variable
  auto   parameters          = dem_parameters;
  auto   physical_properties = dem_parameters.physical_properties;
  double rayleigh_time_step  = 0;

  for (unsigned int i = 0; i < physical_properties.particle_type_number; ++i)
    rayleigh_time_step =
      std::max(M_PI_2 * physical_properties.particle_average_diameter[i] *
                 sqrt(2 * physical_properties.density[i] *
                      (2 + physical_properties.poisson_ratio_particle[i]) *
                      (1 - physical_properties.poisson_ratio_particle[i]) /
                      physical_properties.youngs_modulus_particle[i]) /
                 (0.1631 * physical_properties.poisson_ratio_particle[i] +
                  0.8766),
               rayleigh_time_step);

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
      if (physical_properties.particle_average_diameter.at(i) -
            standard_deviation_multiplier *
              physical_properties.particle_size_std.at(i) <
          0)
        {
          std::cout
            << "Warning: Requested particle size distribution for type: " << i
            << " is not well-defined. Using requested distribution may lead to "
               "changes in sampled particle sizes to avoid negative particle "
               "diameters. You can consider decreasing the standard deviation of "
               "size distribution or increasing the particle diameter instead."
            << std::endl;
        }
    }
}

template void
input_parameter_inspection(const DEMSolverParameters<2> &dem_parameters,
                           const ConditionalOStream &    pcout,
                           const double &standard_deviation_multiplier);

template void
input_parameter_inspection(const DEMSolverParameters<3> &dem_parameters,
                           const ConditionalOStream &    pcout,
                           const double &standard_deviation_multiplier);
