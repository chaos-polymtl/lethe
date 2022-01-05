#include <core/rheological_model.h>

template <int dim>
std::shared_ptr<RheologicalModel<dim>>
RheologicalModel<dim>::model_cast(
  const Parameters::PhysicalProperties &physical_properties)
{
  if (physical_properties.non_newtonian_parameters.model ==
      Parameters::NonNewtonian::Model::powerlaw)
    return std::make_shared<PowerLaw<dim>>(
      physical_properties.non_newtonian_parameters);
  else if (physical_properties.non_newtonian_parameters.model ==
           Parameters::NonNewtonian::Model::carreau)
    return std::make_shared<Carreau<dim>>(
      physical_properties.non_newtonian_parameters);
  else // if (!physical_properties.non_newtonian_flow)
    return std::make_shared<Newtonian<dim>>(
      physical_properties.fluids[0].viscosity);
}

template <int dim>
double 
RheologicalModel<dim>::get_shear_rate_magnitude(const Tensor<2, dim> shear_rate)
{
  double shear_rate_magnitude = 0;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        {
          shear_rate_magnitude += (shear_rate[i][j] * shear_rate[j][i]);
        }
    }
  shear_rate_magnitude = sqrt(0.5 * shear_rate_magnitude);
  return shear_rate_magnitude;
}

template <int dim>
double
Newtonian<dim>::get_viscosity(const double &)
{
  return viscosity;
}

template <int dim>
double
PowerLaw<dim>::get_viscosity(const double &shear_rate_magnitude)
{
  return shear_rate_magnitude > shear_rate_min ?
            K * std::pow(shear_rate_magnitude, n - 1) :
            K * std::pow(shear_rate_min, n - 1);
}

template <int dim>
double
Carreau<dim>::get_viscosity(const double &shear_rate_magnitude)
{
  return viscosity_inf +
          (viscosity_0 - viscosity_inf) *
            std::pow(1.0 + std::pow(shear_rate_magnitude * lambda, a),
                    (n - 1) / a);
}

template class RheologicalModel<2>;
template class RheologicalModel<3>;
template class Newtonian<2>;
template class Newtonian<3>;
template class PowerLaw<2>;
template class PowerLaw<3>;
template class Carreau<2>;
template class Carreau<3>;
