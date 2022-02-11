#include <core/rheological_model.h>

std::shared_ptr<RheologicalModel>
RheologicalModel::model_cast(const Parameters::Fluid &fluid_properties)
{
  if (fluid_properties.rheological_model ==
      Parameters::Fluid::RheologicalModel::newtonian)
    return std::make_shared<Newtonian>(fluid_properties.viscosity);
  else if (fluid_properties.rheological_model ==
           Parameters::Fluid::RheologicalModel::powerlaw)
    return std::make_shared<PowerLaw>(
      fluid_properties.non_newtonian_parameters.powerlaw_parameters.K,
      fluid_properties.non_newtonian_parameters.powerlaw_parameters.n,
      fluid_properties.non_newtonian_parameters.powerlaw_parameters
        .shear_rate_min);
  else // if (physical_properties.fluid[0].non_newtonian_parameters.model ==
       //  Parameters::NonNewtonian::Model::carreau)
    return std::make_shared<Carreau>(
      fluid_properties.non_newtonian_parameters.carreau_parameters.viscosity_0,
      fluid_properties.non_newtonian_parameters.carreau_parameters
        .viscosity_inf,
      fluid_properties.non_newtonian_parameters.carreau_parameters.lambda,
      fluid_properties.non_newtonian_parameters.carreau_parameters.a,
      fluid_properties.non_newtonian_parameters.carreau_parameters.n);
}

double
Newtonian::value(const std::map<field, double> & /*field_values*/)
{
  return viscosity;
}

double
PowerLaw::value(const std::map<field, double> &field_values)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_viscosity(shear_rate_magnitude);
}

void
PowerLaw::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double> &                       property_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    property_vector[i] = calculate_viscosity(shear_rate_magnitude[i]);
}

double
PowerLaw::jacobian(const std::map<field, double> &field_values, const field id)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);
  if (id == field::shear_rate)
    return calculate_derivative(shear_rate_magnitude);
  else
    return 0;
}

void
PowerLaw::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double> &                       jacobian_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  if (id == field::shear_rate)
    for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
      jacobian_vector[i] = calculate_derivative(shear_rate_magnitude[i]);
  else
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
}

double
Carreau::value(const std::map<field, double> &field_values)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_viscosity(shear_rate_magnitude);
}

void
Carreau::vector_value(const std::map<field, std::vector<double>> &field_vectors,
                      std::vector<double> &property_vector)
{
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    {
      property_vector[i] = calculate_viscosity(shear_rate_magnitude[i]);
    }
}

double
Carreau::jacobian(const std::map<field, double> &field_values, const field id)
{
  if (id == field::shear_rate)
    return this->numerical_jacobian(field_values, field::shear_rate);
  else
    return 0;
}

void
Carreau::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double> &                       jacobian_vector)
{
  if (id == field::shear_rate)
    this->vector_numerical_jacobian(field_vectors,
                                    field::shear_rate,
                                    jacobian_vector);
  else
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
}
