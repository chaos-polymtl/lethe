#include <core/rheological_model.h>

template <int dim>
std::shared_ptr<RheologicalModel<dim>>
RheologicalModel<dim>::model_cast(
  const Parameters::PhysicalProperties &physical_properties)
{
  if (!physical_properties.non_newtonian_flow)
    return std::make_shared<Newtonian<dim>>(
      physical_properties.fluids[0].viscosity);
  else if (physical_properties.non_newtonian_parameters.model ==
           Parameters::NonNewtonian::Model::powerlaw)
    return std::make_shared<PowerLaw<dim>>(
      physical_properties.non_newtonian_parameters.powerlaw_parameters.K,
      physical_properties.non_newtonian_parameters.powerlaw_parameters.n,
      physical_properties.non_newtonian_parameters.powerlaw_parameters
        .shear_rate_min);
  else // if (physical_properties.non_newtonian_parameters.model ==
       //  Parameters::NonNewtonian::Model::carreau)
    return std::make_shared<Carreau<dim>>(
      physical_properties.non_newtonian_parameters.carreau_parameters
        .viscosity_0,
      physical_properties.non_newtonian_parameters.carreau_parameters
        .viscosity_inf,
      physical_properties.non_newtonian_parameters.carreau_parameters.lambda,
      physical_properties.non_newtonian_parameters.carreau_parameters.n,
      physical_properties.non_newtonian_parameters.carreau_parameters.a);
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
Newtonian<dim>::value(const std::map<field, double> & /*field_values*/)
{
  return viscosity;
}

template <int dim>
double
PowerLaw<dim>::value(const std::map<field, double> &field_values)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_viscosity(shear_rate_magnitude);
}



template <int dim>
void
PowerLaw<dim>::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double>                        &property_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    property_vector[i] = calculate_viscosity(shear_rate_magnitude[i]);
}

template <int dim>
double
PowerLaw<dim>::jacobian(const std::map<field, double> &field_values,
                        const field                    id)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);
  if (id == field::shear_rate)
    return calculate_derivative(shear_rate_magnitude);
  else
    return 0;
}


template <int dim>
void
PowerLaw<dim>::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  if (id == field::shear_rate)
    for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
      jacobian_vector[i] = calculate_derivative(shear_rate_magnitude[i]);
  else
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
}


template <int dim>
double
Carreau<dim>::value(const std::map<field, double> &field_values)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_viscosity(shear_rate_magnitude);
}

template <int dim>
void
Carreau<dim>::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double>                        &property_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    {
      property_vector[i] = calculate_viscosity(shear_rate_magnitude[i]);
    }
}

template <int dim>
double
Carreau<dim>::jacobian(const std::map<field, double> &field_values,
                       const field                    id)
{
  if (id == field::shear_rate)
    return this->numerical_jacobian(field_values, field::shear_rate);
  else
    return 0;
}


template <int dim>
void
Carreau<dim>::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  if (id == field::shear_rate)
    this->vector_numerical_jacobian(field_vectors,
                                    field::shear_rate,
                                    jacobian_vector);
  else
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
}


template class RheologicalModel<2>;
template class RheologicalModel<3>;
template class Newtonian<2>;
template class Newtonian<3>;
template class PowerLaw<2>;
template class PowerLaw<3>;
template class Carreau<2>;
template class Carreau<3>;
