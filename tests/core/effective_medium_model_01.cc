// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the effective medium models of the electromagnetic properties
 * of a mixture of a fluid (phase 1) and an unresolved dispersed phase
 * (phase 2). The Bruggeman model is checked against reference values. The none
 * model and the field dependencies of the models are also tested. For each
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

  // Values smaller than this threshold are printed as zero, so the output does
  // not depend on the round-off or on the sign of zero
  constexpr double print_threshold = 1e-12;

  using EffectiveMediumModelType =
    Parameters::MaterialInteractions::EffectiveMediumModelType;

  /**
   * @brief Mixture of a fluid and a dispersed phase at a given void fraction
   * and temperature, with its expected effective properties.
   */
  struct MixtureCase
  {
    /// Description of the case, printed in the output
    std::string description;

    /// Material of the fluid (phase 1)
    Parameters::Material fluid;

    /// Material of the dispersed phase (phase 2)
    Parameters::Material dispersed_phase;

    /// Void fraction of the mixture
    double void_fraction;

    /// Temperature of the mixture, only used by temperature-dependent phases
    double temperature;

    /// Expected effective permittivity, which includes the conductivity
    std::complex<double> expected_permittivity;

    /// Expected effective permeability
    std::complex<double> expected_permeability;
  };

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
   * @brief Prints a computed complex value and its expected value, and checks
   * that they agree. The tolerance is relative to the magnitude of the
   * expected value, or absolute if this magnitude is smaller than one.
   *
   * @param[in] computed Computed value.
   * @param[in] expected Expected value.
   */
  void
  check_complex(const std::complex<double> computed,
                const std::complex<double> expected)
  {
    deallog << "  computed : " << format_complex(computed) << std::endl;
    deallog << "  expected : " << format_complex(expected) << std::endl;

    AssertThrow(std::abs(computed - expected) <=
                  tolerance * std::max(std::abs(expected), 1.),
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
   * @brief Tests the Bruggeman model against reference values. Each case is a
   * mixture at a given void fraction and temperature, for which the effective
   * permittivity and permeability are compared with their expected values.
   */
  void
  test_bruggeman_values()
  {
    deallog << "Test 1: Bruggeman model against reference values" << std::endl;

    // Phases with constant properties
    const Parameters::Material lossless_fluid =
      create_material(1., 0., 0., 1., 0.);
    const Parameters::Material lossy_fluid =
      create_material(1., 0.3, 0., 1., 0.);
    const Parameters::Material conductive_dispersed_phase =
      create_material(6., 0.3, 0.1, 2., 0.5);
    const Parameters::Material lossy_dispersed_phase =
      create_material(4., 0.7, 0., 1., 0.);
    const Parameters::Material lossless_dispersed_phase =
      create_material(3., 0., 0., 1., 0.);

    // Phases of the time_harmonic_maxwell_waveguide_bruggeman application test
    // of lethe-fluid-vans
    const Parameters::Material application_test_fluid =
      create_material(1., 0.01, 0.02, 1., 0.01);
    const Parameters::Material application_test_dispersed_phase =
      create_material(0.99687958527608012,
                      0.055119615897750246,
                      0.1,
                      0.99987501701157344,
                      0.035001458134864979);

    // Fluid with a conductivity of 0.001 T. The polynomial coefficients are
    // given from the highest to the lowest degree.
    Parameters::Material temperature_dependent_fluid =
      create_material(1., 0., 0., 1., 0.);
    temperature_dependent_fluid.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;
    temperature_dependent_fluid
      .electric_conductivity_polynomial_coefficients = {0.001, 0.};

    // Dispersed phase with a permittivity of (5 + 0.02 T) + (0.2 + 0.004 T) i,
    // a conductivity of 0.05 + 0.001 T and a permeability of
    // (2 - 0.004 T) + (0.3 + 0.002 T) i
    Parameters::Material temperature_dependent_dispersed_phase =
      create_material(1., 0., 0., 1., 0.);
    temperature_dependent_dispersed_phase.electric_permittivity_model =
      Parameters::Material::ElectricPermittivityModel::polynomial;
    temperature_dependent_dispersed_phase
      .electric_permittivity_real_polynomial_coefficients = {0.02, 5.};
    temperature_dependent_dispersed_phase
      .electric_permittivity_imag_polynomial_coefficients = {0.004, 0.2};
    temperature_dependent_dispersed_phase.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;
    temperature_dependent_dispersed_phase
      .electric_conductivity_polynomial_coefficients = {0.001, 0.05};
    temperature_dependent_dispersed_phase.magnetic_permeability_model =
      Parameters::Material::MagneticPermeabilityModel::polynomial;
    temperature_dependent_dispersed_phase
      .magnetic_permeability_real_polynomial_coefficients = {-0.004, 2.};
    temperature_dependent_dispersed_phase
      .magnetic_permeability_imag_polynomial_coefficients = {0.002, 0.3};

    // The expected values are the admissible roots (positive real part and
    // non-negative imaginary part) of the Bruggeman equation, computed
    // independently in high precision with the complex permittivity
    // eps' + i (eps'' + sigma) of each phase. At the limits of the void
    // fraction, they are the properties of the corresponding phase. With a
    // lossless phase, the exact effective property is real at the
    // corresponding limit, and the computed imaginary part of the root can be
    // slightly negative because of the round-off: these cases check that the
    // root is still accepted. The dispersed phase of the application test was
    // obtained by inverting the Bruggeman equation, so that the mixture is the
    // dissipative medium of the
    // time_harmonic_maxwell_waveguide_robin_dissipative application test of
    // lethe-fluid.
    const std::vector<MixtureCase> mixture_cases = {
      {"Lossless fluid and conductive dispersed phase",
       lossless_fluid,
       conductive_dispersed_phase,
       0.45,
       0.,
       {3.0951335586185476, 0.1481741163000876},
       {1.5009683342308319, 0.22688182043510282}},
      {"Phases of the time_harmonic_maxwell_waveguide_bruggeman application "
       "test",
       application_test_fluid,
       application_test_dispersed_phase,
       0.6,
       0.,
       {1., 0.08},
       {1., 0.02}},
      {"Lossless fluid, limit of the fluid",
       lossless_fluid,
       lossy_dispersed_phase,
       1.,
       0.,
       {1., 0.},
       {1., 0.}},
      {"Lossless fluid, limit of the dispersed phase",
       lossless_fluid,
       lossy_dispersed_phase,
       0.,
       0.,
       {4., 0.7},
       {1., 0.}},
      {"Lossless fluid, clamping of the void fraction above 1",
       lossless_fluid,
       lossy_dispersed_phase,
       1.2,
       0.,
       {1., 0.},
       {1., 0.}},
      {"Lossless fluid, clamping of the void fraction below 0",
       lossless_fluid,
       lossy_dispersed_phase,
       -0.2,
       0.,
       {4., 0.7},
       {1., 0.}},
      {"Lossless dispersed phase, limit of the dispersed phase",
       lossy_fluid,
       lossless_dispersed_phase,
       0.,
       0.,
       {3., 0.},
       {1., 0.}},
      {"Lossless dispersed phase, limit of the fluid",
       lossy_fluid,
       lossless_dispersed_phase,
       1.,
       0.,
       {1., 0.3},
       {1., 0.}},
      {"Temperature-dependent phases, temperature 80",
       temperature_dependent_fluid,
       temperature_dependent_dispersed_phase,
       0.45,
       80.,
       {3.3154167041863749, 0.30986992214162518},
       {1.3528870571729605, 0.21778630219076942}}};

    for (const MixtureCase &mixture_case : mixture_cases)
      {
        const Parameters::PhysicalProperties mixture =
          create_mixture(mixture_case.fluid,
                         mixture_case.dispersed_phase,
                         EffectiveMediumModelType::bruggeman,
                         EffectiveMediumModelType::bruggeman);

        deallog << mixture_case.description << ", void fraction "
                << mixture_case.void_fraction << std::endl;

        deallog << "Permittivity" << std::endl;
        check_complex(
          effective_property(mixture,
                             EffectiveMediumProperty::electric_permittivity,
                             mixture_case.void_fraction,
                             mixture_case.temperature),
          mixture_case.expected_permittivity);

        deallog << "Permeability" << std::endl;
        check_complex(
          effective_property(mixture,
                             EffectiveMediumProperty::magnetic_permeability,
                             mixture_case.void_fraction,
                             mixture_case.temperature),
          mixture_case.expected_permeability);
      }

    deallog << std::endl;
  }

  /**
   * @brief Tests that the vector evaluation of the Bruggeman model gives the
   * same values as its point evaluation.
   */
  void
  test_bruggeman_vector_value()
  {
    deallog << "Test 2: vector evaluation of the Bruggeman model" << std::endl;

    const Parameters::PhysicalProperties mixture =
      create_mixture(create_material(1., 0., 0., 1., 0.),
                     create_material(6., 0.3, 0.1, 2., 0.5),
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
                             void_fractions[i]));
      }

    deallog << std::endl;
  }

  /**
   * @brief Tests the none model, which returns the properties of the fluid.
   */
  void
  test_none_model()
  {
    deallog << "Test 3: none model" << std::endl;

    // The fluid has a permittivity of 2 + 0.1 i, a conductivity of 0.05 and a
    // permeability of 1.5 + 0.2 i
    const Parameters::PhysicalProperties mixture =
      create_mixture(create_material(2., 0.1, 0.05, 1.5, 0.2),
                     create_material(6., 0.3, 0.1, 2., 0.5),
                     EffectiveMediumModelType::none,
                     EffectiveMediumModelType::none);

    deallog << "Permittivity, void fraction " << 0.3 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::electric_permittivity,
                         0.3),
      {2., 0.15});

    deallog << "Permeability, void fraction " << 0.3 << std::endl;
    check_complex(
      effective_property(mixture,
                         EffectiveMediumProperty::magnetic_permeability,
                         0.3),
      {1.5, 0.2});

    deallog << std::endl;
  }

  /**
   * @brief Tests the field dependencies of the effective medium models. The
   * Bruggeman model depends on a field if one of the models of its two phases
   * depends on it, and it also depends on the void fraction. The none model
   * only depends on the fields on which the models of the fluid depend.
   */
  void
  test_dependencies()
  {
    deallog << "Test 4: field dependencies of the models" << std::endl;

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

    // The none model does not use the temperature-dependent dispersed phase
    const Parameters::PhysicalProperties none_mixture =
      create_mixture(create_material(1., 0., 0., 1., 0.),
                     dispersed_phase,
                     EffectiveMediumModelType::none,
                     EffectiveMediumModelType::none);
    const auto none_permittivity_real = EffectiveMediumModel::model_cast_real(
      none_mixture, 0, EffectiveMediumProperty::electric_permittivity);

    deallog << "None permittivity (real part) depends on the void fraction"
            << std::endl;
    check_boolean(none_permittivity_real->depends_on(field::void_fraction),
                  false);

    deallog << "None permittivity (real part) depends on the temperature"
            << std::endl;
    check_boolean(none_permittivity_real->depends_on(field::temperature),
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
      test_bruggeman_values();
      test_bruggeman_vector_value();
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
