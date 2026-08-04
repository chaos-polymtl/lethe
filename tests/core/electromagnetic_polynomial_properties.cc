// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the polynomial temperature-dependent electromagnetic
 *        properties: conductivity, permittivity and permeability.
 */

// C++
#include <cmath>
#include <map>
#include <string>
#include <vector>

// Lethe
#include <core/electric_conductivity_model.h>
#include <core/electric_permittivity_model.h>
#include <core/magnetic_permeability_model.h>

#include <../tests/tests.h>

namespace
{
  constexpr double tolerance = 1e-12;


  template <typename ModelType>
  void
  check_model(const std::string         &name,
              ModelType           &model,
              const std::vector<double> &temperatures,
              const std::vector<double> &expected_values)
  {
    deallog << "Testing " << name << std::endl;

    AssertThrow(temperatures.size() == expected_values.size(),
                ExcMessage("Temperatures and expected values must have the "
                           "same size."));

    for (unsigned int i = 0; i < temperatures.size(); ++i)
      {
        std::map<field, double> field_values;
        field_values[field::temperature] = temperatures[i];

        const double value = model.value(field_values);

        deallog << "  T = " << temperatures[i] << " -> " << value << std::endl;

        AssertThrow(std::abs(value - expected_values[i]) < tolerance,
                    ExcMessage("Unexpected value for " + name));
      }
  }


  void
  test_electric_conductivity()
  {
    deallog << "Testing electric conductivity" << std::endl;

    Parameters::Material material;

    material.electric_conductivity_model =
      Parameters::Material::ElectricConductivityModel::polynomial;

    // Real part:
    // sigma(T) = 3*T^2 - 5
    material.electric_conductivity_polynomial_coefficients = {3.0, 0.0, -5.0};

    const auto model = ElectricConductivityModel::model_cast(material);

    check_model("Electric conductivity",
                *model,
                {0.0, 1.0, -2.0},
                {-5.0, -2.0, 7.0});
  }


  void
  test_electric_permittivity()
  {
    deallog << "Testing electric permittivity" << std::endl;

    Parameters::Material material;

    material.electric_permittivity_model =
      Parameters::Material::ElectricPermittivityModel::polynomial;

    // Real part:
    // epsilon'(T) = 3*T^2 - 5
    material.electric_permittivity_real_polynomial_coefficients = {3.0,
                                                                   0,
                                                                   -5.0};

    // Imaginary part:
    // epsilon''(T) = T^3 + 4*T^2 + T + 1
    material.electric_permittivity_imag_polynomial_coefficients = {1.0,
                                                                   4.0,
                                                                   1.0,
                                                                   1.0};

    const auto real_model =
      ElectricPermittivityModel::model_cast_real(material);

    const auto imag_model =
      ElectricPermittivityModel::model_cast_imag(material);

    check_model("Electric permittivity real",
                *real_model,
                {0.0, 1.0, -2.0},
                {-5.0, -2.0, 7.0});

    check_model("Electric permittivity imaginary",
                *imag_model,
                {0.0, 1.0, -2.0},
                {1.0, 7.0, 7.0});
  }


  void
  test_magnetic_permeability()
  {
    deallog << "Testing magnetic permeability" << std::endl;

    Parameters::Material material;

    material.magnetic_permeability_model =
      Parameters::Material::MagneticPermeabilityModel::polynomial;

    // Real part:
    // mu'(T) = 2
    material.magnetic_permeability_real_polynomial_coefficients = {2.0};

    // Imaginary part:
    // mu''(T) = 2*T^5 - 4*T + 3
    material.magnetic_permeability_imag_polynomial_coefficients = {
      2.0, 0.0, 0.0, 0.0, -4.0, 3.0};

    const auto real_model =
      MagneticPermeabilityModel::model_cast_real(material);

    const auto imag_model =
      MagneticPermeabilityModel::model_cast_imag(material);

    check_model("Magnetic permeability real",
                *real_model,
                {0.0, 1.0, 2.0},
                {2.0, 2.0, 2.0});

    check_model("Magnetic permeability imaginary",
                *imag_model,
                {0.0, 1.0, 2.0},
                {3.0, 1.0, 59.0});
  }


  void
  test()
  {
    deallog << "Begin test" << std::endl;

    test_electric_conductivity();
    test_electric_permittivity();
    test_magnetic_permeability();

    deallog << "OK" << std::endl;
  }

} // namespace


int
main()
{
  try
    {
      initlog();
      test();
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

  return 0;
}
