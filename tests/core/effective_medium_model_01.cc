// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the effective medium models of the electromagnetic properties
 * of a mixture of a fluid (phase 1) and an unresolved dispersed phase
 * (phase 2). The Bruggeman model is checked against reference values computed
 * independently from the closed-form roots of the Bruggeman equation, at the
 * limits of the void fraction (including lossless phases, for which the
 * imaginary part of the computed root is affected by round-off), for the
 * clamping of the void fraction, for the vector evaluation, for the jacobian
 * and for temperature-dependent phase properties. The none model, the field
 * dependencies of the models and the invalid inputs are also tested. For each
 * check, the computed and the expected values are printed.
 */

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <exception>
#include <ios>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Lethe
#include <core/effective_medium_model.h>
#include <core/parameters.h>
#include <core/physical_property_model.h>

// Tests (with common definitions)
#include <../tests/tests.h>

namespace
{
  // Relative tolerance on the comparison with the expected values
  constexpr double tolerance = 1e-12;

  // Relative tolerance on the comparison of the numerical jacobian, which is
  // computed with a first-order forward finite difference, with the analytical
  // derivative
  constexpr double jacobian_tolerance = 1e-5;

  // Values smaller than this threshold are printed as zero, so the output does
  // not depend on the round-off or on the sign of zero
  constexpr double print_threshold = 1e-12;

  using EffectiveMediumModelType =
    Parameters::MaterialInteractions::EffectiveMediumModelType;

  /**
   * @brief Creates a material with constant electromagnetic properties.
   *
   * @param[in] permittivity_real Real part of the relative permittivity.
   * @param[in] permittivity_imag Imaginary part of the relative permittivity.
   * @param[in] conductivity Dimensionless electric conductivity.
   * @param[in] permeability_real Real part of the relative permeability.
   * @param[in] permeability_imag Imaginary part of the relative permeability.
   *
   * @return The material.
   */
  Parameters::Material
  create_material(const double permittivity_real,
                  const double permittivity_imag,
                  const double conductivity,
                  const double permeability_real,
                  const double permeability_imag)
  {
    Parameters::Material material;
    material.electric_permittivity_real = permittivity_real;
    material.electric_permittivity_imag = permittivity_imag;
    material.electric_conductivity      = conductivity;
    material.magnetic_permeability_real = permeability_real;
    material.magnetic_permeability_imag = permeability_imag;
    return material;
  }

  /**
   * @brief Creates the physical properties of a mixture of one fluid and one
   * unresolved dispersed phase, defined by the fluid-dispersed phase material
   * interaction 0.
   *
   * @param[in] fluid Material of the fluid (phase 1).
   * @param[in] dispersed_phase Material of the dispersed phase (phase 2).
   * @param[in] permittivity_model Effective medium model of the permittivity.
   * @param[in] permeability_model Effective medium model of the permeability.
   *
   * @return The physical properties.
   */
  Parameters::PhysicalProperties
  create_mixture(const Parameters::Material    &fluid,
                 const Parameters::Material    &dispersed_phase,
                 const EffectiveMediumModelType permittivity_model,
                 const EffectiveMediumModelType permeability_model)
  {
    Parameters::PhysicalProperties physical_properties;
    physical_properties.number_of_fluids                      = 1;
    physical_properties.fluids                                = {fluid};
    physical_properties.number_of_unresolved_dispersed_phases = 1;
    physical_properties.unresolved_dispersed_phases = {dispersed_phase};

    physical_properties.number_of_material_interactions = 1;
    physical_properties.material_interactions.resize(1);
    Parameters::MaterialInteractions &interaction =
      physical_properties.material_interactions[0];
    interaction.material_interaction_type = Parameters::MaterialInteractions::
      MaterialInteractionsType::fluid_dispersed_phase;
    // Fluid 0 and dispersed phase 0 interact in the material interaction 0
    interaction.fluid_dispersed_phase_interaction_with_material_interaction_id
      .first = {0, 0};
    interaction.fluid_dispersed_phase_interaction_with_material_interaction_id
      .second = 0;
    interaction.electric_permittivity_effective_medium_model =
      permittivity_model;
    interaction.magnetic_permeability_effective_medium_model =
      permeability_model;

    return physical_properties;
  }

  /**
   * @brief Evaluates the complex effective property of the mixture with the
   * models of its real and imaginary parts.
   *
   * @param[in] physical_properties Physical properties of the mixture.
   * @param[in] property Effective property to evaluate.
   * @param[in] void_fraction Void fraction of the mixture.
   * @param[in] temperature Temperature of the mixture. It is only used by the
   * temperature-dependent models of the phases.
   *
   * @return The complex effective property.
   */
  std::complex<double>
  effective_property(const Parameters::PhysicalProperties &physical_properties,
                     const EffectiveMediumProperty         property,
                     const double                          void_fraction,
                     const double                          temperature = 0.)
  {
    std::map<field, double> field_values;
    field_values[field::void_fraction] = void_fraction;
    field_values[field::temperature]   = temperature;

    const auto real_model =
      EffectiveMediumModel::model_cast_real(physical_properties, 0, property);
    const auto imag_model =
      EffectiveMediumModel::model_cast_imag(physical_properties, 0, property);

    return {real_model->value(field_values), imag_model->value(field_values)};
  }

  /**
   * @brief Calculates the derivative of the complex effective property of the
   * mixture with respect to a field, with the jacobians of the models of its
   * real and imaginary parts.
   *
   * @param[in] physical_properties Physical properties of the mixture.
   * @param[in] property Effective property to differentiate.
   * @param[in] void_fraction Void fraction of the mixture.
   * @param[in] temperature Temperature of the mixture.
   * @param[in] field_id Field with respect to which the derivative is
   * calculated.
   *
   * @return The derivative of the complex effective property.
   */
  std::complex<double>
  effective_property_jacobian(
    const Parameters::PhysicalProperties &physical_properties,
    const EffectiveMediumProperty         property,
    const double                          void_fraction,
    const double                          temperature,
    const field                           field_id)
  {
    std::map<field, double> field_values;
    field_values[field::void_fraction] = void_fraction;
    field_values[field::temperature]   = temperature;

    const auto real_model =
      EffectiveMediumModel::model_cast_real(physical_properties, 0, property);
    const auto imag_model =
      EffectiveMediumModel::model_cast_imag(physical_properties, 0, property);

    return {real_model->jacobian(field_values, field_id),
            imag_model->jacobian(field_values, field_id)};
  }

  /**
   * @brief Calculates the derivative of the Bruggeman effective property
   * alpha_eff with respect to a variable x, by implicit differentiation of the
   * Bruggeman equation F = epsilon g_1 + (1 - epsilon) g_2 = 0, with
   * g_k = (alpha_k - alpha_eff) / (alpha_k + 2 alpha_eff):
   * d alpha_eff / dx = - (dF/d epsilon d epsilon/dx
   *                       + dF/d alpha_1 d alpha_1/dx
   *                       + dF/d alpha_2 d alpha_2/dx) / (dF/d alpha_eff),
   * where
   * dF/d epsilon   = g_1 - g_2,
   * dF/d alpha_1   = 3 epsilon alpha_eff / (alpha_1 + 2 alpha_eff)^2,
   * dF/d alpha_2   = 3 (1 - epsilon) alpha_eff / (alpha_2 + 2 alpha_eff)^2,
   * dF/d alpha_eff = -3 (epsilon alpha_1 / (alpha_1 + 2 alpha_eff)^2
   *                      + (1 - epsilon) alpha_2 / (alpha_2 + 2 alpha_eff)^2).
   *
   * @param[in] void_fraction Void fraction epsilon of the mixture.
   * @param[in] alpha_1 Property of phase 1.
   * @param[in] alpha_2 Property of phase 2.
   * @param[in] alpha_eff Effective property of the mixture.
   * @param[in] d_void_fraction Derivative of the void fraction with respect to
   * x.
   * @param[in] d_alpha_1 Derivative of the property of phase 1 with respect to
   * x.
   * @param[in] d_alpha_2 Derivative of the property of phase 2 with respect to
   * x.
   *
   * @return The derivative of the effective property with respect to x.
   */
  std::complex<double>
  bruggeman_derivative(const double               void_fraction,
                       const std::complex<double> alpha_1,
                       const std::complex<double> alpha_2,
                       const std::complex<double> alpha_eff,
                       const double               d_void_fraction,
                       const std::complex<double> d_alpha_1,
                       const std::complex<double> d_alpha_2)
  {
    const std::complex<double> denominator_1 = alpha_1 + 2. * alpha_eff;
    const std::complex<double> denominator_2 = alpha_2 + 2. * alpha_eff;
    const std::complex<double> g_1 = (alpha_1 - alpha_eff) / denominator_1;
    const std::complex<double> g_2 = (alpha_2 - alpha_eff) / denominator_2;

    const std::complex<double> dF_dvoid_fraction = g_1 - g_2;
    const std::complex<double> dF_dalpha_1 =
      3. * void_fraction * alpha_eff / (denominator_1 * denominator_1);
    const std::complex<double> dF_dalpha_2 =
      3. * (1. - void_fraction) * alpha_eff / (denominator_2 * denominator_2);
    const std::complex<double> dF_dalpha_eff =
      -3. * (void_fraction * alpha_1 / (denominator_1 * denominator_1) +
             (1. - void_fraction) * alpha_2 / (denominator_2 * denominator_2));

    return -(dF_dvoid_fraction * d_void_fraction + dF_dalpha_1 * d_alpha_1 +
             dF_dalpha_2 * d_alpha_2) /
           dF_dalpha_eff;
  }

  /**
   * @brief Returns the value to print, with the values below the print
   * threshold set to zero.
   *
   * @param[in] value Value to print.
   *
   * @return The printed value.
   */
  double
  printed(const double value)
  {
    return (std::abs(value) < print_threshold) ? 0. : value;
  }

  /**
   * @brief Formats a complex value as "a + b i", with the number format of
   * deallog (6 significant digits, trailing zeros shown).
   *
   * @param[in] value Complex value to format.
   *
   * @return The formatted value.
   */
  std::string
  format_complex(const std::complex<double> value)
  {
    const double real_part = printed(std::real(value));
    const double imag_part = printed(std::imag(value));

    std::ostringstream stream;
    stream << std::showpoint << real_part << (imag_part < 0. ? " - " : " + ")
           << std::abs(imag_part) << " i";
    return stream.str();
  }

  /**
   * @brief Prints the constant electromagnetic properties of a phase.
   *
   * @param[in] name Name of the phase.
   * @param[in] material Material of the phase.
   */
  void
  print_phase(const std::string &name, const Parameters::Material &material)
  {
    deallog << name << ": permittivity "
            << format_complex({material.electric_permittivity_real,
                               material.electric_permittivity_imag})
            << ", conductivity " << material.electric_conductivity
            << ", permeability "
            << format_complex({material.magnetic_permeability_real,
                               material.magnetic_permeability_imag})
            << std::endl;
  }

  /**
   * @brief Prints a computed complex value and its expected value, and checks
   * that they agree. The tolerance is relative to the magnitude of the
   * expected value, or absolute if this magnitude is smaller than one.
   *
   * @param[in] computed Computed value.
   * @param[in] expected Expected value.
   * @param[in] relative_tolerance Tolerance of the comparison.
   */
  void
  check_complex(const std::complex<double> computed,
                const std::complex<double> expected,
                const double               relative_tolerance)
  {
    deallog << "  computed : " << format_complex(computed) << std::endl;
    deallog << "  expected : " << format_complex(expected) << std::endl;

    AssertThrow(std::abs(computed - expected) <=
                  relative_tolerance * std::max(std::abs(expected), 1.),
                ExcMessage("The computed value differs from the expected "
                           "value."));
  }

  /**
   * @brief Prints a computed boolean and its expected value, and checks that
   * they are equal.
   *
   * @param[in] computed Computed boolean.
   * @param[in] expected Expected boolean.
   */
  void
  check_boolean(const bool computed, const bool expected)
  {
    deallog << "  computed : " << (computed ? "true" : "false") << std::endl;
    deallog << "  expected : " << (expected ? "true" : "false") << std::endl;

    AssertThrow(computed == expected,
                ExcMessage("The computed boolean differs from the expected "
                           "boolean."));
  }

  /**
   * @brief Evaluates the real part of the effective permittivity and returns
   * true if an exception is thrown, either when the model is created or when
   * it is evaluated.
   *
   * @param[in] physical_properties Physical properties of the mixture.
   * @param[in] field_values Values of the fields.
   *
   * @return True if an exception is thrown.
   */
  bool
  evaluation_throws(const Parameters::PhysicalProperties &physical_properties,
                    const std::map<field, double>        &field_values)
  {
    try
      {
        const auto model = EffectiveMediumModel::model_cast_real(
          physical_properties,
          0,
          EffectiveMediumProperty::electric_permittivity);
        model->value(field_values);
      }
    catch (const std::exception &)
      {
        return true;
      }
    return false;
  }

  /**
   * @brief Tests the Bruggeman model against reference values. The reference
   * values are the admissible roots of the Bruggeman equation, computed
   * independently.
   */
  void
  test_bruggeman_reference_values()
  {
    deallog << "Test 1: Bruggeman model against reference values" << std::endl;

    const Parameters::Material fluid = create_material(1., 0., 0., 1., 0.);
    const Parameters::Material dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    print_phase("Fluid", fluid);
    print_phase("Dispersed phase", dispersed_phase);

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    deallog << "Permittivity, void fraction " << 0.45 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         0.45),
      {3.0951335586185476, 0.1481741163000876},
      tolerance);

    deallog << "Permeability, void fraction " << 0.3 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::magnetic_permeability,
                         0.3),
      {1.659430424772245, 0.31132070834320985},
      tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests the Bruggeman model with the properties of the
   * time_harmonic_maxwell_waveguide_bruggeman application test of
   * lethe-fluid-vans. The properties of the dispersed phase were obtained by
   * inverting the Bruggeman equation, so the mixture has the properties of the
   * dissipative medium of the time_harmonic_maxwell_waveguide_robin_dissipative
   * application test of lethe-fluid.
   */
  void
  test_bruggeman_application_test_mixture()
  {
    deallog << "Test 2: Bruggeman model with the properties of the "
               "time_harmonic_maxwell_waveguide_bruggeman application test"
            << std::endl;

    const Parameters::Material fluid =
      create_material(1., 0.01, 0.02, 1., 0.01);
    const Parameters::Material dispersed_phase =
      create_material(0.99687958527608012,
                      0.055119615897750246,
                      0.1,
                      0.99987501701157344,
                      0.035001458134864979);
    print_phase("Fluid", fluid);
    print_phase("Dispersed phase", dispersed_phase);

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    deallog << "Permittivity, void fraction " << 0.6 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         0.6),
      {1., 0.08},
      tolerance);

    deallog << "Permeability, void fraction " << 0.6 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::magnetic_permeability,
                         0.6),
      {1., 0.02},
      tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests the Bruggeman model at the limits of the void fraction and
   * the clamping of the void fraction. With a lossless phase, the exact
   * effective property is real at the corresponding limit, and the computed
   * imaginary part of the root can be slightly negative because of the
   * round-off. These cases test that the root is still accepted.
   */
  void
  test_bruggeman_limits()
  {
    deallog << "Test 3: Bruggeman model at the limits of the void fraction"
            << std::endl;

    // Lossless fluid: the effective property is real for a void fraction of 1
    const Parameters::Material lossless_fluid =
      create_material(1., 0., 0., 1., 0.);
    const Parameters::Material lossy_dispersed_phase =
      create_material(4., 0.7, 0., 1., 0.);
    print_phase("Fluid", lossless_fluid);
    print_phase("Dispersed phase", lossy_dispersed_phase);

    const Parameters::PhysicalProperties lossless_fluid_mixture =
      create_mixture(lossless_fluid,
                     lossy_dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    const std::vector<double> void_fractions = {1., 0., 1.2, -0.2};
    const std::vector<std::complex<double>> expected_values = {{1., 0.},
                                                               {4., 0.7},
                                                               {1., 0.},
                                                               {4., 0.7}};
    for (unsigned int i = 0; i < void_fractions.size(); ++i)
      {
        deallog << "Permittivity, void fraction " << void_fractions[i]
                << std::endl;
        check_complex(
          effective_property(lossless_fluid_mixture,
                             EffectiveMediumProperty::electric_permittivity,
                             void_fractions[i]),
          expected_values[i],
          tolerance);
      }

    // Lossless dispersed phase: the effective property is real for a void
    // fraction of 0
    const Parameters::Material lossy_fluid =
      create_material(1., 0.3, 0., 1., 0.);
    const Parameters::Material lossless_dispersed_phase =
      create_material(3., 0., 0., 1., 0.);
    print_phase("Fluid", lossy_fluid);
    print_phase("Dispersed phase", lossless_dispersed_phase);

    const Parameters::PhysicalProperties lossless_dispersed_phase_mixture =
      create_mixture(lossy_fluid,
                     lossless_dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    deallog << "Permittivity, void fraction " << 0. << std::endl;
    check_complex(
      effective_property(lossless_dispersed_phase_mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         0.),
      {3., 0.},
      tolerance);

    deallog << "Permittivity, void fraction " << 1. << std::endl;
    check_complex(
      effective_property(lossless_dispersed_phase_mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         1.),
      {1., 0.3},
      tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests that the vector evaluation of the Bruggeman model gives the
   * same values as its point evaluation.
   */
  void
  test_bruggeman_vector_value()
  {
    deallog << "Test 4: vector evaluation of the Bruggeman model" << std::endl;

    const Parameters::Material fluid = create_material(1., 0., 0., 1., 0.);
    const Parameters::Material dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    print_phase("Fluid", fluid);
    print_phase("Dispersed phase", dispersed_phase);

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    const std::vector<double> void_fractions = {0., 0.3, 0.45, 0.8, 1.};
    const unsigned int        n_points       = void_fractions.size();

    std::map<field, std::vector<double>> field_vectors;
    field_vectors[field::void_fraction] = void_fractions;

    const auto real_model = EffectiveMediumModel::model_cast_real(
      mixture, 0, EffectiveMediumProperty::electric_permittivity);
    const auto imag_model = EffectiveMediumModel::model_cast_imag(
      mixture, 0, EffectiveMediumProperty::electric_permittivity);

    std::vector<double> real_values(n_points);
    std::vector<double> imag_values(n_points);
    real_model->vector_value(field_vectors, real_values);
    imag_model->vector_value(field_vectors, imag_values);

    for (unsigned int i = 0; i < n_points; ++i)
      {
        deallog << "Permittivity, void fraction " << void_fractions[i]
                << std::endl;
        check_complex(
          {real_values[i], imag_values[i]},
          effective_property(mixture,
                             EffectiveMediumProperty::electric_permittivity,
                             void_fractions[i]),
          tolerance);
      }

    deallog << std::endl;
  }

  /**
   * @brief Tests the jacobian of the Bruggeman model with respect to the void
   * fraction against the analytical derivative obtained by implicit
   * differentiation of the Bruggeman equation, and its jacobian with respect to
   * the temperature for phases with constant properties.
   */
  void
  test_bruggeman_jacobian()
  {
    deallog << "Test 5: jacobian of the Bruggeman model" << std::endl;

    const Parameters::Material fluid = create_material(1., 0., 0., 1., 0.);
    const Parameters::Material dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    print_phase("Fluid", fluid);
    print_phase("Dispersed phase", dispersed_phase);

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    // Permittivities of the phases and reference effective permittivity of
    // test 1
    const double               void_fraction = 0.45;
    const std::complex<double> alpha_1       = {1., 0.};
    const std::complex<double> alpha_2       = {6., 0.4};
    const std::complex<double> alpha_eff     = {3.0951335586185476,
                                                0.1481741163000876};

    deallog << "Derivative of the permittivity with respect to the void "
               "fraction, void fraction "
            << void_fraction << std::endl;
    check_complex(effective_property_jacobian(
                    mixture,
                    EffectiveMediumProperty::electric_permittivity,
                    void_fraction,
                    0.,
                    field::void_fraction),
                  bruggeman_derivative(
                    void_fraction, alpha_1, alpha_2, alpha_eff, 1., 0., 0.),
                  jacobian_tolerance);

    // The properties of the phases do not depend on the temperature
    deallog << "Derivative of the permittivity with respect to the "
               "temperature, void fraction "
            << void_fraction << std::endl;
    check_complex(effective_property_jacobian(
                    mixture,
                    EffectiveMediumProperty::electric_permittivity,
                    void_fraction,
                    0.,
                    field::temperature),
                  {0., 0.},
                  jacobian_tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests the Bruggeman model with phase properties that depend on the
   * temperature. The temperature is used by the models of both phases. The
   * expected values are the admissible roots of the Bruggeman equation,
   * computed independently with the properties of the phases at the given
   * temperature, and the analytical derivative with respect to the
   * temperature.
   */
  void
  test_bruggeman_temperature_dependent_phases()
  {
    deallog << "Test 6: Bruggeman model with temperature-dependent phase "
               "properties"
            << std::endl;

    // Fluid with a conductivity of 0.001 T. The polynomial coefficients are
    // given from the highest to the lowest degree.
    Parameters::Material fluid = create_material(1., 0., 0., 1., 0.);
    fluid.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;
    fluid.electric_conductivity_polynomial_coefficients = {0.001, 0.};

    // Dispersed phase with a permittivity of (5 + 0.02 T) + (0.2 + 0.004 T) i,
    // a conductivity of 0.05 + 0.001 T and a permeability of
    // (2 - 0.004 T) + (0.3 + 0.002 T) i
    Parameters::Material dispersed_phase = create_material(1., 0., 0., 1., 0.);
    dispersed_phase.electric_permittivity_model =
      Parameters::Material::ElectricPermittivityModel::polynomial;
    dispersed_phase.electric_permittivity_real_polynomial_coefficients = {0.02,
                                                                          5.};
    dispersed_phase.electric_permittivity_imag_polynomial_coefficients = {0.004,
                                                                          0.2};
    dispersed_phase.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;
    dispersed_phase.electric_conductivity_polynomial_coefficients = {0.001,
                                                                     0.05};
    dispersed_phase.magnetic_permeability_model =
      Parameters::Material::MagneticPermeabilityModel::polynomial;
    dispersed_phase.magnetic_permeability_real_polynomial_coefficients = {
      -0.004, 2.};
    dispersed_phase.magnetic_permeability_imag_polynomial_coefficients = {0.002,
                                                                          0.3};

    deallog << "Fluid: permittivity 1, conductivity 0.001 T, permeability 1"
            << std::endl;
    deallog << "Dispersed phase: permittivity (5 + 0.02 T) + (0.2 + 0.004 T) "
               "i, conductivity 0.05 + 0.001 T, permeability (2 - 0.004 T) + "
               "(0.3 + 0.002 T) i"
            << std::endl;

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    const double                            void_fraction = 0.45;
    const std::vector<double>               temperatures  = {20., 80.};
    const std::vector<std::complex<double>> expected_permittivities = {
      {2.8708712356656223, 0.14838611386677864},
      {3.3154167041863749, 0.30986992214162518}};
    const std::vector<std::complex<double>> expected_permeabilities = {
      {1.4609719811603226, 0.15612159802650313},
      {1.3528870571729605, 0.21778630219076942}};

    // Point evaluation
    for (unsigned int i = 0; i < temperatures.size(); ++i)
      {
        deallog << "Permittivity, temperature " << temperatures[i]
                << ", void fraction " << void_fraction << std::endl;
        check_complex(
          effective_property(mixture,
                             EffectiveMediumProperty::electric_permittivity,
                             void_fraction,
                             temperatures[i]),
          expected_permittivities[i],
          tolerance);

        deallog << "Permeability, temperature " << temperatures[i]
                << ", void fraction " << void_fraction << std::endl;
        check_complex(
          effective_property(mixture,
                             EffectiveMediumProperty::magnetic_permeability,
                             void_fraction,
                             temperatures[i]),
          expected_permeabilities[i],
          tolerance);
      }

    // Vector evaluation
    const unsigned int                   n_points = temperatures.size();
    std::map<field, std::vector<double>> field_vectors;
    field_vectors[field::void_fraction] =
      std::vector<double>(n_points, void_fraction);
    field_vectors[field::temperature] = temperatures;

    const auto real_model = EffectiveMediumModel::model_cast_real(
      mixture, 0, EffectiveMediumProperty::electric_permittivity);
    const auto imag_model = EffectiveMediumModel::model_cast_imag(
      mixture, 0, EffectiveMediumProperty::electric_permittivity);

    std::vector<double> real_values(n_points);
    std::vector<double> imag_values(n_points);
    real_model->vector_value(field_vectors, real_values);
    imag_model->vector_value(field_vectors, imag_values);

    for (unsigned int i = 0; i < n_points; ++i)
      {
        deallog << "Permittivity (vector evaluation), temperature "
                << temperatures[i] << ", void fraction " << void_fraction
                << std::endl;
        check_complex({real_values[i], imag_values[i]},
                      expected_permittivities[i],
                      tolerance);
      }

    // Derivative of the permittivity with respect to the temperature. The
    // permittivities of the phases (with their conductivity) and their
    // derivatives with respect to the temperature are:
    // alpha_1 = 1 + 0.001 T i, d alpha_1 / dT = 0.001 i,
    // alpha_2 = (5 + 0.02 T) + (0.25 + 0.005 T) i,
    // d alpha_2 / dT = 0.02 + 0.005 i.
    const double               temperature = temperatures[0];
    const std::complex<double> alpha_1     = {1., 0.001 * temperature};
    const std::complex<double> alpha_2     = {5. + 0.02 * temperature,
                                              0.25 + 0.005 * temperature};
    const std::complex<double> d_alpha_1   = {0., 0.001};
    const std::complex<double> d_alpha_2   = {0.02, 0.005};

    deallog << "Derivative of the permittivity with respect to the "
               "temperature, temperature "
            << temperature << ", void fraction " << void_fraction << std::endl;
    check_complex(effective_property_jacobian(
                    mixture,
                    EffectiveMediumProperty::electric_permittivity,
                    void_fraction,
                    temperature,
                    field::temperature),
                  bruggeman_derivative(void_fraction,
                                       alpha_1,
                                       alpha_2,
                                       expected_permittivities[0],
                                       0.,
                                       d_alpha_1,
                                       d_alpha_2),
                  jacobian_tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests the none model, which returns the properties of the fluid.
   */
  void
  test_none_model()
  {
    deallog << "Test 7: none model" << std::endl;

    const Parameters::Material fluid = create_material(2., 0.1, 0.05, 1.5, 0.2);
    const Parameters::Material dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    print_phase("Fluid", fluid);
    print_phase("Dispersed phase", dispersed_phase);

    const Parameters::PhysicalProperties mixture =
      create_mixture(fluid,
                     dispersed_phase,
                     EffectiveMediumModelType::none,
                     EffectiveMediumModelType::none);

    deallog << "Permittivity, void fraction " << 0.3 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         0.3),
      {2., 0.15},
      tolerance);

    deallog << "Permeability, void fraction " << 0.3 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::magnetic_permeability,
                         0.3),
      {1.5, 0.2},
      tolerance);

    deallog << std::endl;
  }

  /**
   * @brief Tests the field dependencies of the effective medium models. A
   * mixture depends on a field if one of the models of its phases depends on
   * it, and the Bruggeman model also depends on the void fraction.
   */
  void
  test_dependencies()
  {
    deallog << "Test 8: field dependencies of the models" << std::endl;

    // Dispersed phase with a temperature-dependent conductivity
    Parameters::Material dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    dispersed_phase.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;
    dispersed_phase.electric_conductivity_polynomial_coefficients = {0.01, 0.1};

    const Parameters::PhysicalProperties bruggeman_mixture =
      create_mixture(create_material(1., 0., 0., 1., 0.),
                     dispersed_phase,
                     EffectiveMediumModelType::bruggeman,
                     EffectiveMediumModelType::bruggeman);

    const auto permittivity_real = EffectiveMediumModel::model_cast_real(
      bruggeman_mixture, 0, EffectiveMediumProperty::electric_permittivity);
    const auto permittivity_imag = EffectiveMediumModel::model_cast_imag(
      bruggeman_mixture, 0, EffectiveMediumProperty::electric_permittivity);
    const auto permeability_real = EffectiveMediumModel::model_cast_real(
      bruggeman_mixture, 0, EffectiveMediumProperty::magnetic_permeability);

    deallog << "Bruggeman permittivity (real part) depends on the void "
               "fraction"
            << std::endl;
    check_boolean(permittivity_real->depends_on(field::void_fraction), true);

    deallog << "Bruggeman permittivity (real part) depends on the temperature"
            << std::endl;
    check_boolean(permittivity_real->depends_on(field::temperature), true);

    deallog << "Bruggeman permittivity (imaginary part) depends on the "
               "temperature"
            << std::endl;
    check_boolean(permittivity_imag->depends_on(field::temperature), true);

    deallog << "Bruggeman permeability (real part) depends on the void "
               "fraction"
            << std::endl;
    check_boolean(permeability_real->depends_on(field::void_fraction), true);

    deallog << "Bruggeman permeability (real part) depends on the temperature"
            << std::endl;
    check_boolean(permeability_real->depends_on(field::temperature), false);

    const Parameters::PhysicalProperties none_mixture =
      create_mixture(create_material(1., 0., 0., 1., 0.),
                     create_material(6., 0.3, 0.1, 2., 0.5),
                     EffectiveMediumModelType::none,
                     EffectiveMediumModelType::none);
    deallog << "None permittivity (real part) depends on the void fraction"
            << std::endl;
    check_boolean(EffectiveMediumModel::model_cast_real(
                    none_mixture,
                    0,
                    EffectiveMediumProperty::electric_permittivity)
                    ->depends_on(field::void_fraction),
                  false);

    deallog << std::endl;
  }
} // namespace

int
main()
{
  try
    {
      initlog();
      test_bruggeman_reference_values();
      test_bruggeman_application_test_mixture();
      test_bruggeman_limits();
      test_bruggeman_vector_value();
      test_bruggeman_jacobian();
      test_bruggeman_temperature_dependent_phases();
      test_none_model();
      test_dependencies();
      deallog << "OK" << std::endl;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
