// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/rheological_model.h>

std::shared_ptr<RheologicalModel>
RheologicalModel::model_cast(const Parameters::Material &material_properties)
{
  if (material_properties.rheological_model ==
      Parameters::Material::RheologicalModel::newtonian)
    return std::make_shared<Newtonian>(material_properties.kinematic_viscosity);
  else if (material_properties.rheological_model ==
           Parameters::Material::RheologicalModel::powerlaw)
    return std::make_shared<PowerLaw>(
      material_properties.non_newtonian_parameters.powerlaw_parameters.K,
      material_properties.non_newtonian_parameters.powerlaw_parameters.n,
      material_properties.non_newtonian_parameters.powerlaw_parameters
        .shear_rate_min);
  else if (material_properties.rheological_model ==
           Parameters::Material::RheologicalModel::carreau)
    return std::make_shared<Carreau>(
      material_properties.non_newtonian_parameters.carreau_parameters
        .kinematic_viscosity_0,
      material_properties.non_newtonian_parameters.carreau_parameters
        .kinematic_viscosity_inf,
      material_properties.non_newtonian_parameters.carreau_parameters.lambda,
      material_properties.non_newtonian_parameters.carreau_parameters.a,
      material_properties.non_newtonian_parameters.carreau_parameters.n);

  else // (fluid_properties.rheological_model ==
       //    Parameters::Fluid::RheologicalModel::phase_change)
    return std::make_shared<PhaseChangeRheology>(
      material_properties.phase_change_parameters);
}

double
Newtonian::value(const std::map<field, double> & /*field_values*/)
{
  return kinematic_viscosity;
}

double
PowerLaw::value(const std::map<field, double> &field_values)
{
  Assert(field_values.find(field::shear_rate) != field_values.end(),
         PhysicialPropertyModelFieldUndefined("PowerLaw", "shear_rate"));
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_kinematic_viscosity(shear_rate_magnitude);
}

void
PowerLaw::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double>                        &property_vector)
{
  Assert(field_vectors.find(field::shear_rate) != field_vectors.end(),
         PhysicialPropertyModelFieldUndefined("PowerLaw", "shear_rate"));
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    property_vector[i] = calculate_kinematic_viscosity(shear_rate_magnitude[i]);
}

double
PowerLaw::jacobian(const std::map<field, double> &field_values, const field id)
{
  const double shear_rate_magnitude = field_values.at(field::shear_rate);
  if (id == field::shear_rate)
    {
      Assert(field_values.find(field::shear_rate) != field_values.end(),
             PhysicialPropertyModelFieldUndefined("PowerLaw", "shear_rate"));
      return calculate_derivative(shear_rate_magnitude);
    }
  else
    return 0;
}

void
PowerLaw::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  Assert(field_vectors.find(field::shear_rate) != field_vectors.end(),
         PhysicialPropertyModelFieldUndefined("PowerLaw", "shear_rate"));
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
  Assert(field_values.find(field::shear_rate) != field_values.end(),
         PhysicialPropertyModelFieldUndefined("Carreau", "shear_rate"));
  const double shear_rate_magnitude = field_values.at(field::shear_rate);

  return calculate_kinematic_viscosity(shear_rate_magnitude);
}

void
Carreau::vector_value(const std::map<field, std::vector<double>> &field_vectors,
                      std::vector<double> &property_vector)
{
  Assert(field_vectors.find(field::shear_rate) != field_vectors.end(),
         PhysicialPropertyModelFieldUndefined("Carreau", "shear_rate"));
  const auto shear_rate_magnitude = field_vectors.at(field::shear_rate);

  for (unsigned int i = 0; i < shear_rate_magnitude.size(); ++i)
    {
      property_vector[i] =
        calculate_kinematic_viscosity(shear_rate_magnitude[i]);
    }
}

double
Carreau::jacobian(const std::map<field, double> &field_values, const field id)
{
  if (id == field::shear_rate)
    {
      Assert(field_values.find(field::shear_rate) != field_values.end(),
             PhysicialPropertyModelFieldUndefined("Carreau", "shear_rate"));
      return this->numerical_jacobian(field_values, field::shear_rate);
    }
  else
    return 0;
}

void
Carreau::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  if (id == field::shear_rate)
    {
      Assert(field_vectors.find(field::shear_rate) != field_vectors.end(),
             PhysicialPropertyModelFieldUndefined("Carreau", "shear_rate"));
      this->vector_numerical_jacobian(field_vectors,
                                      field::shear_rate,
                                      jacobian_vector);
    }
  else
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
}



double
PhaseChangeRheology::value(const std::map<field, double> &field_values)
{
  Assert(field_values.find(field::temperature) != field_values.end(),
         PhysicialPropertyModelFieldUndefined("PhaseChangeRheology",
                                              "temperature"));
  const double temperature = field_values.at(field::temperature);

  return kinematic_viscosity(temperature);
}

void
PhaseChangeRheology::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double>                        &property_vector)
{
  Assert(field_vectors.find(field::temperature) != field_vectors.end(),
         PhysicialPropertyModelFieldUndefined("PhaseChangeRheology",
                                              "temperature"));
  const std::vector<double> &temperature_vec =
    field_vectors.at(field::temperature);

  for (unsigned int i = 0; i < temperature_vec.size(); ++i)
    {
      property_vector[i] = kinematic_viscosity(temperature_vec[i]);
    }
}

double
PhaseChangeRheology::jacobian(const std::map<field, double> &field_values,
                              const field                    id)
{
  if (id == field::temperature)
    {
      Assert(field_values.find(field::temperature) != field_values.end(),
             PhysicialPropertyModelFieldUndefined("PhaseChangeRheology",
                                                  "temperature"));
      return this->numerical_jacobian(field_values, field::temperature);
    }
  else
    return 0;
}

void
PhaseChangeRheology::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  Assert(field_vectors.find(field::temperature) != field_vectors.end(),
         PhysicialPropertyModelFieldUndefined("PhaseChangeRheology",
                                              "temperature"));
  vector_numerical_jacobian(field_vectors, id, jacobian_vector);
}

/**
 * @brief Calculates the kinematic viscosity used in PSPG and SUPG stabilizations.
 * @param[in] field_values Value of the various fields on which the property may
 * depend.
 */
double
PhaseChangeRheology::get_kinematic_viscosity_for_stabilization(
  const std::map<field, double> & /*field_values*/)
{
  return param.kinematic_viscosity_l;
}

/**
 * @brief Calculates the vector values of the kinematic viscosity used in PSPG and SUPG stabilizations.
 * @param[in] field_vectors Value of the fields on which the property may
 * depend on.
 * @param[out] property_vector Vector of computed viscosities.
 */
void
PhaseChangeRheology::get_kinematic_viscosity_for_stabilization_vector(
  const std::map<field, std::vector<double>> & /*field_vectors*/,
  std::vector<double> &property_vector)
{
  std::fill(property_vector.begin(),
            property_vector.end(),
            param.kinematic_viscosity_l);
}

/**
 * @brief Calculates the dynamic viscosity used in PSPG and SUPG stabilizations.
 * @param[in] p_density_ref The density of the fluid at the reference state
 * @param[in] field_values Value of the various fields on which the property may
 * depend.
 */
double
PhaseChangeRheology::get_dynamic_viscosity_for_stabilization(
  const double &p_density_ref,
  const std::map<field, double> & /*field_values*/)
{
  return param.kinematic_viscosity_l * p_density_ref;
}

/**
 * @brief Calculates the vector values of the dynamic viscosity used in PSPG and SUPG stabilizations.
 * @param[in] p_density_ref The density of the fluid at the reference state
 * @param[in] field_vectors Value of the fields on which the property may
 * depend on.
 * @param[out] property_vector Vector of computed viscosities.
 */
void
PhaseChangeRheology::get_dynamic_viscosity_for_stabilization_vector(
  const double &p_density_ref,
  const std::map<field, std::vector<double>> & /*field_vectors*/,
  std::vector<double> &property_vector)
{
  std::fill(property_vector.begin(),
            property_vector.end(),
            param.kinematic_viscosity_l * p_density_ref);
}
