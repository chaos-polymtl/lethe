// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/effective_medium_model.h>
#include <core/electric_conductivity_model.h>
#include <core/electric_permittivity_model.h>
#include <core/magnetic_permeability_model.h>

#include <algorithm>
#include <cmath>

EffectiveMediumModel::EffectiveMediumModel(
  const Parameters::PhysicalProperties &physical_properties,
  const unsigned int                    material_interaction_id,
  const EffectiveMediumProperty         property,
  const Component                       component)
  : property(property)
  , component(component)
{
  AssertThrow(material_interaction_id <
                physical_properties.number_of_material_interactions,
              ExcMessage("The material interaction of the effective medium "
                         "model does not exist."));

  const Parameters::MaterialInteractions &interaction =
    physical_properties.material_interactions[material_interaction_id];

  // The material interaction gives the materials of the two phases of the
  // mixture: the fluid and the unresolved particles for a fluid-particles
  // interaction, the two fluids for a fluid-fluid interaction.
  if (interaction.material_interaction_type ==
      Parameters::MaterialInteractions::MaterialInteractionsType::
        fluid_particles)
    {
      const auto &[fluid_id, particles_id] =
        interaction.fluid_particles_interaction_with_material_interaction_id
          .first;
      AssertThrow(fluid_id < physical_properties.number_of_fluids &&
                    particles_id <
                      physical_properties.number_of_unresolved_particles,
                  ExcMessage("The fluid or the particles of the "
                             "fluid-particles interaction do not exist."));

      phase_1_models =
        create_phase_models(physical_properties.fluids[fluid_id], property);
      phase_2_models = create_phase_models(
        physical_properties.unresolved_particles[particles_id], property);
    }
  else if (interaction.material_interaction_type ==
           Parameters::MaterialInteractions::MaterialInteractionsType::
             fluid_fluid)
    {
      const auto &[fluid_0_id, fluid_1_id] =
        interaction.fluid_fluid_interaction_with_material_interaction_id.first;
      AssertThrow(fluid_0_id < physical_properties.number_of_fluids &&
                    fluid_1_id < physical_properties.number_of_fluids,
                  ExcMessage("The fluids of the fluid-fluid interaction do "
                             "not exist."));

      phase_1_models =
        create_phase_models(physical_properties.fluids[fluid_0_id], property);
      phase_2_models =
        create_phase_models(physical_properties.fluids[fluid_1_id], property);
    }
  else
    AssertThrow(false,
                ExcMessage("Effective medium models are only available for "
                           "fluid-particles and fluid-fluid interactions."));

  // The mixture property depends on a field if at least one model of either
  // phase depends on it (union of the phase model dependencies).
  for (auto &[field_id, mixture_depends_on_field] : this->model_depends_on)
    {
      for (const auto &model : phase_1_models)
        mixture_depends_on_field =
          mixture_depends_on_field || model->depends_on(field_id);
      for (const auto &model : phase_2_models)
        mixture_depends_on_field =
          mixture_depends_on_field || model->depends_on(field_id);
    }
}

std::shared_ptr<EffectiveMediumModel>
EffectiveMediumModel::model_cast_real(
  const Parameters::PhysicalProperties &physical_properties,
  const unsigned int                    material_interaction_id,
  const EffectiveMediumProperty         property)
{
  const Parameters::MaterialInteractions &interaction =
    physical_properties.material_interactions.at(material_interaction_id);

  // Assign the effective medium model type based on the property being
  // calculated
  Parameters::MaterialInteractions::EffectiveMediumModelType model_type;
  if (property == EffectiveMediumProperty::electric_permittivity)
    model_type = interaction.electric_effective_medium_model;
  else if (property == EffectiveMediumProperty::magnetic_permeability)
    model_type = interaction.magnetic_effective_medium_model;
  else
    AssertThrow(false,
                ExcMessage("Invalid effective medium property. The choices are "
                           "<electric_permittivity|magnetic_permeability>"));

  if (model_type ==
      Parameters::MaterialInteractions::EffectiveMediumModelType::none)
    return std::make_shared<NoEffectiveMediumModel>(physical_properties,
                                                    material_interaction_id,
                                                    property,
                                                    Component::real);
  else if (model_type == Parameters::MaterialInteractions::
                           EffectiveMediumModelType::bruggeman)
    return std::make_shared<BruggemanEffectiveMediumModel>(
      physical_properties, material_interaction_id, property, Component::real);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid effective medium model. The choices are <none|bruggeman>"));
}

std::shared_ptr<EffectiveMediumModel>
EffectiveMediumModel::model_cast_imag(
  const Parameters::PhysicalProperties &physical_properties,
  const unsigned int                    material_interaction_id,
  const EffectiveMediumProperty         property)
{
  const Parameters::MaterialInteractions &interaction =
    physical_properties.material_interactions.at(material_interaction_id);

  // Assign the effective medium model type based on the property being
  // calculated
  Parameters::MaterialInteractions::EffectiveMediumModelType model_type;
  if (property == EffectiveMediumProperty::electric_permittivity)
    model_type = interaction.electric_effective_medium_model;
  else if (property == EffectiveMediumProperty::magnetic_permeability)
    model_type = interaction.magnetic_effective_medium_model;
  else
    AssertThrow(false,
                ExcMessage("Invalid effective medium property. The choices are "
                           "<electric_permittivity|magnetic_permeability>"));

  if (model_type ==
      Parameters::MaterialInteractions::EffectiveMediumModelType::none)
    return std::make_shared<NoEffectiveMediumModel>(physical_properties,
                                                    material_interaction_id,
                                                    property,
                                                    Component::imag);
  else if (model_type == Parameters::MaterialInteractions::
                           EffectiveMediumModelType::bruggeman)
    return std::make_shared<BruggemanEffectiveMediumModel>(
      physical_properties, material_interaction_id, property, Component::imag);
  else
    AssertThrow(
      false,
      ExcMessage(
        "Invalid effective medium model. The choices are <none|bruggeman>"));
}

double
EffectiveMediumModel::value(const std::map<field, double> &field_values)
{
  const std::complex<double> phase_1_property =
    phase_value(phase_1_models, field_values);
  const std::complex<double> phase_2_property =
    phase_value(phase_2_models, field_values);

  return get_component(
    this->apply_mixing_law(field_values, phase_1_property, phase_2_property));
}

void
EffectiveMediumModel::vector_value(
  const std::map<field, std::vector<double>> &field_vectors,
  std::vector<double>                        &property_vector)
{
  const unsigned int n_points = property_vector.size();

  // Local containers: the model is shared between threads, so it must not
  // store intermediate results.
  std::vector<std::complex<double>> phase_1_property(n_points);
  std::vector<std::complex<double>> phase_2_property(n_points);
  std::vector<std::complex<double>> mixture_property(n_points);

  phase_vector_value(phase_1_models, field_vectors, phase_1_property);
  phase_vector_value(phase_2_models, field_vectors, phase_2_property);

  this->vector_apply_mixing_law(field_vectors,
                                phase_1_property,
                                phase_2_property,
                                mixture_property);

  for (unsigned int i = 0; i < n_points; ++i)
    property_vector[i] = get_component(mixture_property[i]);
}

double
EffectiveMediumModel::jacobian(const std::map<field, double> &field_values,
                               const field                    id)
{
  return this->numerical_jacobian(field_values, id);
}

void
EffectiveMediumModel::vector_jacobian(
  const std::map<field, std::vector<double>> &field_vectors,
  const field                                 id,
  std::vector<double>                        &jacobian_vector)
{
  this->vector_numerical_jacobian(field_vectors, id, jacobian_vector);
}

std::vector<std::shared_ptr<PhysicalPropertyModel>>
EffectiveMediumModel::create_phase_models(
  const Parameters::Material   &material_properties,
  const EffectiveMediumProperty property)
{
  std::vector<std::shared_ptr<PhysicalPropertyModel>> phase_models;

  switch (property)
    {
      case EffectiveMediumProperty::electric_permittivity:
        // Models of eps', eps'' and sigma
        phase_models.push_back(
          ElectricPermittivityModel::model_cast_real(material_properties));
        phase_models.push_back(
          ElectricPermittivityModel::model_cast_imag(material_properties));
        phase_models.push_back(
          ElectricConductivityModel::model_cast(material_properties));
        break;
      case EffectiveMediumProperty::magnetic_permeability:
        // Models of mu' and mu''
        phase_models.push_back(
          MagneticPermeabilityModel::model_cast_real(material_properties));
        phase_models.push_back(
          MagneticPermeabilityModel::model_cast_imag(material_properties));
        break;
    }

  return phase_models;
}

std::complex<double>
EffectiveMediumModel::phase_value(
  const std::vector<std::shared_ptr<PhysicalPropertyModel>> &phase_models,
  const std::map<field, double>                             &field_values) const
{
  switch (property)
    {
      case EffectiveMediumProperty::electric_permittivity:
        // Complex permittivity: eps' + i (eps'' + sigma)
        return {phase_models[0]->value(field_values),
                phase_models[1]->value(field_values) +
                  phase_models[2]->value(field_values)};
      case EffectiveMediumProperty::magnetic_permeability:
        // Complex permeability: mu' + i mu''
        return {phase_models[0]->value(field_values),
                phase_models[1]->value(field_values)};
    }

  AssertThrow(false, ExcMessage("Invalid effective medium property."));
  return {0.0, 0.0};
}

void
EffectiveMediumModel::phase_vector_value(
  const std::vector<std::shared_ptr<PhysicalPropertyModel>> &phase_models,
  const std::map<field, std::vector<double>>                &field_vectors,
  std::vector<std::complex<double>> &phase_property) const
{
  const unsigned int  n_points = phase_property.size();
  std::vector<double> real_part(n_points);
  std::vector<double> imaginary_part(n_points);

  switch (property)
    {
      case EffectiveMediumProperty::electric_permittivity:
        {
          // Complex permittivity: eps' + i (eps'' + sigma)
          std::vector<double> conductivity(n_points);
          phase_models[0]->vector_value(field_vectors, real_part);
          phase_models[1]->vector_value(field_vectors, imaginary_part);
          phase_models[2]->vector_value(field_vectors, conductivity);
          for (unsigned int i = 0; i < n_points; ++i)
            imaginary_part[i] += conductivity[i];
          break;
        }
      case EffectiveMediumProperty::magnetic_permeability:
        // Complex permeability: mu' + i mu''
        phase_models[0]->vector_value(field_vectors, real_part);
        phase_models[1]->vector_value(field_vectors, imaginary_part);
        break;
    }

  for (unsigned int i = 0; i < n_points; ++i)
    phase_property[i] = {real_part[i], imaginary_part[i]};
}

double
EffectiveMediumModel::get_component(
  const std::complex<double> mixture_property) const
{
  return (component == Component::real) ? std::real(mixture_property) :
                                          std::imag(mixture_property);
}

BruggemanEffectiveMediumModel::BruggemanEffectiveMediumModel(
  const Parameters::PhysicalProperties &physical_properties,
  const unsigned int                    material_interaction_id,
  const EffectiveMediumProperty         property,
  const Component                       component)
  : EffectiveMediumModel(physical_properties,
                         material_interaction_id,
                         property,
                         component)
{
  // The Bruggeman mixing law depends on the void fraction
  this->model_depends_on[field::void_fraction] = true;
}

std::complex<double>
BruggemanEffectiveMediumModel::apply_mixing_law(
  const std::map<field, double> &field_values,
  const std::complex<double>     phase_1_property,
  const std::complex<double>     phase_2_property) const
{
  AssertThrow(field_values.contains(field::void_fraction),
              PhysicialPropertyModelFieldUndefined(
                "BruggemanEffectiveMediumModel", "void_fraction"));

  return bruggeman(field_values.at(field::void_fraction),
                   phase_1_property,
                   phase_2_property);
}

void
BruggemanEffectiveMediumModel::vector_apply_mixing_law(
  const std::map<field, std::vector<double>> &field_vectors,
  const std::vector<std::complex<double>>    &phase_1_property,
  const std::vector<std::complex<double>>    &phase_2_property,
  std::vector<std::complex<double>>          &mixture_property) const
{
  AssertThrow(field_vectors.contains(field::void_fraction),
              PhysicialPropertyModelFieldUndefined(
                "BruggemanEffectiveMediumModel", "void_fraction"));

  const std::vector<double> &void_fraction =
    field_vectors.at(field::void_fraction);
  AssertDimension(void_fraction.size(), mixture_property.size());
  AssertDimension(phase_1_property.size(), mixture_property.size());
  AssertDimension(phase_2_property.size(), mixture_property.size());

  for (unsigned int i = 0; i < mixture_property.size(); ++i)
    mixture_property[i] =
      bruggeman(void_fraction[i], phase_1_property[i], phase_2_property[i]);
}

std::complex<double>
BruggemanEffectiveMediumModel::bruggeman(
  const double               void_fraction,
  const std::complex<double> phase_1_property,
  const std::complex<double> phase_2_property)
{
  AssertThrow(std::real(phase_1_property) > 0.0 &&
                std::imag(phase_1_property) >= 0.0,
              InvalidPhysicalPropertyComplex(phase_1_property));
  AssertThrow(std::real(phase_2_property) > 0.0 &&
                std::imag(phase_2_property) >= 0.0,
              InvalidPhysicalPropertyComplex(phase_2_property));

  // We clamp the void fraction to be between 0 and 1, since the Bruggeman model
  // is only valid for void fractions in this range and slightly out of range
  // values can occur from the projected void fraction due to numerical errors.
  const double phase_1_fraction = std::clamp(void_fraction, 0., 1.);
  const double phase_2_fraction = 1.0 - phase_1_fraction;
  const std::complex<double> factor =
    (3.0 * phase_2_fraction - 1.0) * phase_2_property +
    (3.0 * phase_1_fraction - 1.0) * phase_1_property;
  const std::complex<double> discriminant_root =
    std::sqrt(factor * factor + 8.0 * phase_2_property * phase_1_property);

  // The std::sqrt function gives one of the branches of the complex plane. To
  // ensure that the effective property is in the upper-half right quadrant of
  // the complex plane, we compute both roots and select the proper one.
  const std::complex<double> root_1 = 0.25 * (factor + discriminant_root);
  const std::complex<double> root_2 = 0.25 * (factor - discriminant_root);

  if (std::real(root_1) > 0.0 && std::imag(root_1) >= 0.0)
    return root_1;
  if (std::real(root_2) > 0.0 && std::imag(root_2) >= 0.0)
    return root_2;

  AssertThrow(
    false,
    ExcMessage(
      "The computed effective property is not in the upper-half right quadrant of the complex plane."));
  return {0.0, 0.0};
}
