#include <dem/input_parameter_inspection.h>

using namespace dealii;

template <int dim>
void
input_parameter_inspection(const DEMSolverParameters<dim> &dem_parameters)
{
  // Getting the input parameters as local variable
  auto   parameters         = dem_parameters;
  double rayleigh_time_step = 0;

  for (unsigned int i = 0;
       i < dem_parameters.physical_properties.particle_type_number;
       ++i)
    rayleigh_time_step = std::max(
      M_PI_2 * parameters.physical_properties.particle_average_diameter[i] *
        sqrt(2 * parameters.physical_properties.density[i] *
             (2 + parameters.physical_properties.poisson_ratio_particle[i]) *
             (1 - parameters.physical_properties.poisson_ratio_particle[i]) /
             parameters.physical_properties.youngs_modulus_particle[i]) /
        (0.1631 * parameters.physical_properties.poisson_ratio_particle[i] +
         0.8766),
      rayleigh_time_step);

  const double time_step_rayleigh_ratio =
    parameters.simulation_control.dt / rayleigh_time_step;
  std::cout << "DEM time-step is " << time_step_rayleigh_ratio * 100
            << "% of Rayleigh time step" << std::endl;

  if (time_step_rayleigh_ratio > 0.15)
    {
      std::cout << "Warning: It is recommended to decrease the time-step"
                << std::endl;
    }
  else if (time_step_rayleigh_ratio < 0.01)
    {
      std::cout << "Warning: It is recommended to increase the time-step"
                << std::endl;
    }
}

template void
input_parameter_inspection(const DEMSolverParameters<2> &dem_parameters);

template void
input_parameter_inspection(const DEMSolverParameters<3> &dem_parameters);
