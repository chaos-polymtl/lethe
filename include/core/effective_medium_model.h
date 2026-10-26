// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_effective_medium_model_h
#define lethe_effective_medium_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

#include <deal.II/base/exceptions.h>

#include <complex>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

DeclException1(
  InvalidPhysicalPropertyComplex,
  std::complex<double>,
  << "Invalid physical property value: " << arg1
  << ". The physical property values must be in the upper-half right quadrant of the complex plane.");

/**
 * @brief Physical property of a mixture that can be calculated with an
 * effective medium model.
 */
enum class EffectiveMediumProperty : std::int8_t
{
  electric_permittivity,
  magnetic_permeability
};

/**
 * @brief Effective medium model for the calculation of a physical property of a mixture
 * of two phases. The mixture and its two phases are defined by a material
 * interaction. The property of each phase is calculated with the physical
 * property models of its material. The mixing is done with complex values to
 * allow for the calculation of complex properties, such as the complex
 * permittivity or permeability. In the case of a real property, the imaginary
 * part is simply set to zero. The field on which the mixing law depends are the
 * ones of the fields on which the phase models depend and the ones required by
 * the mixing law.
 */
class EffectiveMediumModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Flag for the complex number part of the mixture property that is returned by the model. This is required because a PhysicalPropertyModel can only return a single double value, but mixing laws for complex properties require the calculation of both the real and imaginary parts in the same model.
   */
  enum class Component : std::int8_t
  {
    real,
    imag
  };

  /**
   * @brief Constructor. Gets the materials of the two phases from the material
   * interaction, creates the models of their property and sets the fields on
   * which these models depend.
   *
   * @param[in] physical_properties Physical properties parameters.
   * @param[in] material_interaction_id Index of the material interaction that
   * defines the mixture.
   * @param[in] property Mixture property calculated by the model.
   * @param[in] component Part of the mixture property returned by the model.
   */
  EffectiveMediumModel(
    const Parameters::PhysicalProperties &physical_properties,
    const unsigned int                    material_interaction_id,
    const EffectiveMediumProperty         property,
    const Component                       component);

  /**
   * @brief Creates the model of the real part of a mixture property, with the
   * mixing law given by the material interaction.
   *
   * @param[in] physical_properties Physical properties parameters.
   * @param[in] material_interaction_id Index of the material interaction that
   * defines the mixture.
   * @param[in] property Mixture property calculated by the model.
   *
   * @return A shared pointer to the created effective medium model.
   */
  static std::shared_ptr<EffectiveMediumModel>
  model_cast_real(const Parameters::PhysicalProperties &physical_properties,
                  const unsigned int                    material_interaction_id,
                  const EffectiveMediumProperty         property);

  /**
   * @brief Creates the model of the imaginary part of a mixture property,
   * with the mixing law given by the material interaction.
   *
   * @param[in] physical_properties Physical properties parameters.
   * @param[in] material_interaction_id Index of the material interaction that
   * defines the mixture.
   * @param[in] property Mixture property calculated by the model.
   *
   * @return A shared pointer to the created effective medium model.
   */
  static std::shared_ptr<EffectiveMediumModel>
  model_cast_imag(const Parameters::PhysicalProperties &physical_properties,
                  const unsigned int                    material_interaction_id,
                  const EffectiveMediumProperty         property);

  /**
   * @brief Calculates the selected part of the mixture property.
   *
   * @param[in] field_values Values of the fields on which the property
   * depends.
   *
   * @return Real or imaginary part of the mixture property.
   */
  double
  value(const std::map<field, double> &field_values) override;

  /**
   * @brief Calculates the selected part of the mixture property at multiple
   * points.
   *
   * @param[in] field_vectors Values of the fields on which the property
   * depends at each point.
   * @param[in,out] property_vector Real or imaginary part of the mixture
   * property at each point. Its size sets the number of points.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override;

  /**
   * @brief Calculates the derivative of the selected part of the mixture
   * property with respect to a field using a finite difference.
   *
   * @param[in] field_values Values of the fields on which the property
   * depends.
   * @param[in] id Field with respect to which the derivative is calculated.
   *
   * @return Derivative of the property with respect to the field.
   */
  double
  jacobian(const std::map<field, double> &field_values,
           const field                    id) override;

  /**
   * @brief Calculates the derivative of the selected part of the mixture
   * property with respect to a field at multiple points using a finite
   * difference.
   *
   * @param[in] field_vectors Values of the fields on which the property
   * depends at each point.
   * @param[in] id Field with respect to which the derivative is calculated.
   * @param[in,out] jacobian_vector Derivative of the property at each point.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override;

  /**
   * @brief Mixing law of the effective medium model.
   *
   * @param[in] field_values Values of the fields on which the mixing law
   * depends.
   * @param[in] phase_1_property Property of phase 1.
   * @param[in] phase_2_property Property of phase 2.
   *
   * @return Property of the mixture.
   */
  virtual std::complex<double>
  apply_mixing_law(const std::map<field, double> &field_values,
                   const std::complex<double>     phase_1_property,
                   const std::complex<double>     phase_2_property) const = 0;

  /**
   * @brief Mixing law of the effective medium model at multiple points.
   *
   * @param[in] field_vectors Values of the fields on which the mixing law
   * depends at each point.
   * @param[in] phase_1_property Property of phase 1 at each point.
   * @param[in] phase_2_property Property of phase 2 at each point.
   * @param[in,out] mixture_property Property of the mixture at each point.
   */
  virtual void
  vector_apply_mixing_law(
    const std::map<field, std::vector<double>> &field_vectors,
    const std::vector<std::complex<double>>    &phase_1_property,
    const std::vector<std::complex<double>>    &phase_2_property,
    std::vector<std::complex<double>>          &mixture_property) const = 0;

private:
  /**
   * @brief Creates the models for the two phases needed to calculate a property of a material.
   *
   * @param[in] material_properties Parameters of the material.
   * @param[in] property Property described by the models.
   *
   * @return Models of the property.
   */
  static std::vector<std::shared_ptr<PhysicalPropertyModel>>
  create_phase_models(const Parameters::Material   &material_properties,
                      const EffectiveMediumProperty property);

  /**
   * @brief Calculates the property of a phase from its models, according to
   * the mixture property.
   *
   * @param[in] phase_models Models of the property of the phase.
   * @param[in] field_values Values of the fields.
   *
   * @return Property of the phase.
   */
  std::complex<double>
  phase_value(
    const std::vector<std::shared_ptr<PhysicalPropertyModel>> &phase_models,
    const std::map<field, double> &field_values) const;

  /**
   * @brief Calculates the property of a phase at multiple points from its
   * models, according to the mixture property.
   *
   * @param[in] phase_models Models of the property of the phase.
   * @param[in] field_vectors Values of the fields at each point.
   * @param[in,out] phase_property Property of the phase at each point. Its
   * size sets the number of points.
   */
  void
  phase_vector_value(
    const std::vector<std::shared_ptr<PhysicalPropertyModel>> &phase_models,
    const std::map<field, std::vector<double>>                &field_vectors,
    std::vector<std::complex<double>> &phase_property) const;

  /**
   * @brief Returns the selected complex component of a mixture property.
   *
   * @param[in] mixture_property Property of the mixture.
   *
   * @return Real or imaginary part of the mixture property.
   */
  double
  get_component(const std::complex<double> mixture_property) const;

  /**
   * @brief Returns the effective medium model chosen in a material interaction
   * for a mixture property.
   *
   * @param[in] interaction Parameters of the material interaction.
   * @param[in] property Mixture property calculated by the model.
   *
   * @return Type of the effective medium model of the property.
   */
  static Parameters::MaterialInteractions::EffectiveMediumModelType
  get_model_type(const Parameters::MaterialInteractions &interaction,
                 const EffectiveMediumProperty           property);

  /// Mixture property calculated by the model
  const EffectiveMediumProperty property;

  /// Part of the mixture property returned by the model
  const Component component;

  /// Models of the property of phase 1
  std::vector<std::shared_ptr<PhysicalPropertyModel>> phase_1_models;

  /// Models of the property of phase 2
  std::vector<std::shared_ptr<PhysicalPropertyModel>> phase_2_models;
};

/**
 * @brief Mixture without effective medium model: the mixture property is the
 * property of phase 1. It is used when one of the phases as no impact on the
 * mixture property.
 */
class NoEffectiveMediumModel final : public EffectiveMediumModel
{
public:
  using EffectiveMediumModel::EffectiveMediumModel;

  /**
   * @brief Returns the property of phase 1.
   *
   * @param[in] field_values Values of the fields (unused).
   * @param[in] phase_1_property Property of phase 1.
   * @param[in] phase_2_property Property of phase 2 (unused).
   *
   * @return Property of phase 1.
   */
  std::complex<double>
  apply_mixing_law(
    [[maybe_unused]] const std::map<field, double> &field_values,
    const std::complex<double>                      phase_1_property,
    [[maybe_unused]] const std::complex<double> phase_2_property) const override
  {
    return phase_1_property;
  }

  /**
   * @brief Returns the property of phase 1 at multiple points.
   *
   * @param[in] field_vectors Values of the fields at each point (unused).
   * @param[in] phase_1_property Property of phase 1 at each point.
   * @param[in] phase_2_property Property of phase 2 at each point (unused).
   * @param[in,out] mixture_property Property of phase 1 at each point.
   */
  void
  vector_apply_mixing_law(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    const std::vector<std::complex<double>>                  &phase_1_property,
    [[maybe_unused]] const std::vector<std::complex<double>> &phase_2_property,
    std::vector<std::complex<double>> &mixture_property) const override
  {
    mixture_property = phase_1_property;
  }
};

/**
 * @brief Two-phase Bruggeman effective medium approximation for spherical
 * inclusions. With the void fraction \f$\epsilon\f$ (volume fraction of phase
 * 1), the mixture property \f$\alpha_\mathrm{eff}\f$ is the root of
 * \f$\epsilon\frac{\alpha_1 - \alpha_\mathrm{eff}}{\alpha_1 +
 * 2\alpha_\mathrm{eff}} + (1-\epsilon)\frac{\alpha_2 -
 * \alpha_\mathrm{eff}}{\alpha_2 + 2\alpha_\mathrm{eff}} = 0\f$, a quadratic
 * equation. The root with a positive real part and a non-negative imaginary
 * part is kept.
 */
class BruggemanEffectiveMediumModel final : public EffectiveMediumModel
{
public:
  /**
   * @brief Constructor. Adds the void fraction to the fields on which the
   * mixture property depends.
   *
   * @param[in] physical_properties Physical properties parameters.
   * @param[in] material_interaction_id Index of the material interaction that
   * defines the mixture.
   * @param[in] property Mixture property calculated by the model.
   * @param[in] component Part of the mixture property returned by the model.
   */
  BruggemanEffectiveMediumModel(
    const Parameters::PhysicalProperties &physical_properties,
    const unsigned int                    material_interaction_id,
    const EffectiveMediumProperty         property,
    const Component                       component);

  /**
   * @brief Calculates the mixture property with the Bruggeman effective medium
   * approximation.
   *
   * @param[in] field_values Values of the fields. It must contain the void
   * fraction.
   * @param[in] phase_1_property Property of phase 1.
   * @param[in] phase_2_property Property of phase 2.
   *
   * @return Property of the mixture.
   */
  std::complex<double>
  apply_mixing_law(const std::map<field, double> &field_values,
                   const std::complex<double>     phase_1_property,
                   const std::complex<double> phase_2_property) const override;

  /**
   * @brief Calculates the mixture property with the Bruggeman effective medium
   * approximation at multiple points.
   *
   * @param[in] field_vectors Values of the fields at each point. It must
   * contain the void fraction.
   * @param[in] phase_1_property Property of phase 1 at each point.
   * @param[in] phase_2_property Property of phase 2 at each point.
   * @param[in,out] mixture_property Property of the mixture at each point.
   */
  void
  vector_apply_mixing_law(
    const std::map<field, std::vector<double>> &field_vectors,
    const std::vector<std::complex<double>>    &phase_1_property,
    const std::vector<std::complex<double>>    &phase_2_property,
    std::vector<std::complex<double>> &mixture_property) const override;

private:
  /**
   * @brief Bruggeman formula for a given void fraction.
   *
   * @param[in] void_fraction Void fraction of the mixture (volume fraction of
   * phase 1).
   * @param[in] phase_1_property Property of phase 1.
   * @param[in] phase_2_property Property of phase 2.
   *
   * @return Property of the mixture.
   */
  static std::complex<double>
  bruggeman(const double               void_fraction,
            const std::complex<double> phase_1_property,
            const std::complex<double> phase_2_property);
};

#endif
