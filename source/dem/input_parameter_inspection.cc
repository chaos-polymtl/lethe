#include <dem/input_parameter_inspection.h>

using namespace dealii;

template <int dim>
void
input_parameter_inspection(const DEMSolverParameters<dim> &dem_parameters)
{
  // Getting the input parameters as local variable
  auto parameters = dem_parameters;

  const double rayleigh_time_step =
    M_PI_2 * parameters.physical_properties.diameter *
    sqrt(2 * parameters.physical_properties.density *
         (2 + parameters.physical_properties.poisson_ratio_particle) *
         (1 - parameters.physical_properties.poisson_ratio_particle) /
         parameters.physical_properties.youngs_modulus_particle) /
    (0.1631 * parameters.physical_properties.poisson_ratio_particle + 0.8766);

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
