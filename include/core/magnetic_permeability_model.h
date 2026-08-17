// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_magnetic_permeability_model_h
#define lethe_magnetic_permeability_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * magnetic permeability on each quadrature point.
 */
class MagneticPermeabilityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a MagneticPermeabilityModel object by casting it to
   * the proper child class for the real part of the permeability.
   *
   * @param[in] material_properties Parameters for a material
   */
  static std::shared_ptr<MagneticPermeabilityModel>
  model_cast_real(const Parameters::Material &material_properties);

  /**
   * @brief Instantiates and returns a pointer to a MagneticPermeabilityModel object by casting it to
   * the proper child class for the imaginary part of the permeability.
   *
   * @param[in] material_properties Parameters for a material
   */
  static std::shared_ptr<MagneticPermeabilityModel>
  model_cast_imag(const Parameters::Material &material_properties);
};


/**
 * @brief Constant magnetic permeability.
 */
class ConstantMagneticPermeability : public MagneticPermeabilityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantMagneticPermeability(const double p_magnetic_permeability)
    : magnetic_permeability(p_magnetic_permeability)
  {}

  /**
   * @brief Calculates the value the magnetic permeability
   * @param[in] field_values Values of the various fields on which the magnetic
   * permeability depends.
   * @return Value of the magnetic permeability calculated with the field_values.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return magnetic_permeability;
  };

  /**
   * @brief Calculates, in a vector, values of magnetic permeabilities
   * @param[in] field_vectors Vector of properties on which the magnetic
   * permeabilities depend
   * @param[in,out] property_vector Values of the magnetic permeabilities
   */
  void
  vector_value(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double> &property_vector) override
  {
    property_vector.assign(property_vector.size(), magnetic_permeability);
  }

  /**
   * @brief Calculates the jacobian (the partial derivative) of the magnetic permeability with respect to a field
   * @param[in] field_values Value of the various fields on which the property
   * may depend.
   * @param[in] id Identifier of the field with respect to which the jacobian
   * should be calculated
   * @return Value of the partial derivative of the magnetic permeability with respect to the specified field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief Calculates, in a vector, the derivatives of the magnetic permeability with respect to a field
   * @param[in] field_vectors Vector for the values of the fields used to
   * evaluate the property
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be calculated
   * @param[in,out] jacobian Vector of values of the derivatives of the magnetic
   * permeability with respect to the specified field
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    [[maybe_unused]] const field id,
    std::vector<double>         &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double magnetic_permeability;
};

/**
 * @brief Polynomial magnetic permeability.
 */
class PolynomialMagneticPermeability : public MagneticPermeabilityModel
{
public:
  /**
   * @brief Default constructor
   */
  PolynomialMagneticPermeability(
    const std::vector<double> &p_magnetic_permeability_polynomial_coefficients)
    : magnetic_permeability_polynomial_coefficients(
        p_magnetic_permeability_polynomial_coefficients)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief Calculates the value the magnetic permeability
   * @param[in] field_values Values of the various fields on which the magnetic
   * permeability depends.
   * @return Value of the magnetic permeability calculated with the field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    AssertThrow(field_values.contains(field::temperature),
                PhysicialPropertyModelFieldUndefined(
                  "PolynomialMagneticPermeability", "temperature"));
    const double temperature = field_values.at(field::temperature);

    // Compute the polynomial power using Horner's method for efficiency, i.e.
    // P(x) = (...((a_n * x + a_{n-1}) * x + a_{n-2}) * x + ... + a_1) * x + a_0
    double magnetic_permeability = 0.0;
    for (const auto &coefficient :
         magnetic_permeability_polynomial_coefficients)
      {
        magnetic_permeability =
          coefficient + temperature * magnetic_permeability;
      }

    return magnetic_permeability;
  }

  /**
   * @brief Calculates, in a vector, values of magnetic permeabilities
   * @param[in] field_vectors Vector of properties on which the magnetic
   * permeabilities depend
   * @param[in,out] property_vector Values of the magnetic permeabilities
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    AssertThrow(field_vectors.find(field::temperature) != field_vectors.end(),
                PhysicialPropertyModelFieldUndefined(
                  "PolynomialMagneticPermeability", "temperature"));

    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      {
        // Compute the polynomial power using Horner's method for efficiency,
        // i.e. P(x) = (...((a_n * x + a_{n-1}) * x + a_{n-2}) * x + ... + a_1)
        // * x + a_0
        double magnetic_permeability = 0.0;
        for (const auto &coefficient :
             magnetic_permeability_polynomial_coefficients)
          {
            magnetic_permeability =
              coefficient + temperature[i] * magnetic_permeability;
          }
        property_vector[i] = magnetic_permeability;
      }
  }

  /**
   * @brief Calculates the jacobian (the partial derivative) of the magnetic permeability with respect to a field
   * @param[in] field_values Value of the various fields on which the property
   * may depend.
   * @param[in] id Identifier of the field with respect to which the jacobian
   * should be calculated
   * @return Value of the partial derivative of the magnetic permeability with respect to the specified field.
   *
   * @remark The jacobian is only implemented with respect to the temperature.
   */
  double
  jacobian(const std::map<field, double> &field_values, const field id) override
  {
    // The derivative can only be taken with respect to the temperature field,
    // since the polynomial is a function of temperature.
    if (id == field::temperature)
      {
        AssertThrow(field_values.find(field::temperature) != field_values.end(),
                    PhysicialPropertyModelFieldUndefined(
                      "PolynomialMagneticPermeability", "temperature"));

        // Use a Horner accumulator to evaluate the derivative directly.
        unsigned int polynomial_order =
          magnetic_permeability_polynomial_coefficients.size() - 1;
        double magnetic_permeability_derivative =
          polynomial_order * magnetic_permeability_polynomial_coefficients[0];

        double temperature = field_values.at(field::temperature);

        for (unsigned int i = 1; i < polynomial_order; ++i)
          {
            magnetic_permeability_derivative =
              magnetic_permeability_derivative * temperature +
              (polynomial_order - i) *
                magnetic_permeability_polynomial_coefficients[i];
          }

        return magnetic_permeability_derivative;
      }
    else
      return 0;
  }

  /**
   * @brief Calculates, in a vector, the derivatives of the magnetic permeability with respect to a field
   * @param[in] field_vectors Vector for the values of the fields used to
   * evaluate the property
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be calculated
   * @param[in,out] jacobian Vector of values of the derivatives of the magnetic
   * permeability with respect to the specified field
   *
   * @remark The vector_jacobian is only implemented with respect to the
   * temperature.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    if (id == field::temperature)
      {
        AssertThrow(
          field_vectors.find(field::temperature) != field_vectors.end(),
          PhysicialPropertyModelFieldUndefined("PolynomialMagneticPermeability",
                                               "temperature"));

        const std::vector<double> &temperature =
          field_vectors.at(field::temperature);
        for (unsigned int i = 0; i < jacobian_vector.size(); ++i)
          {
            // Use a Horner accumulator to evaluate the derivative directly.
            unsigned int polynomial_order =
              magnetic_permeability_polynomial_coefficients.size() - 1;
            double magnetic_permeability_derivative =
              polynomial_order *
              magnetic_permeability_polynomial_coefficients[0];
            for (unsigned int j = 1; j < polynomial_order; ++j)
              {
                magnetic_permeability_derivative =
                  magnetic_permeability_derivative * temperature[i] +
                  (polynomial_order - j) *
                    magnetic_permeability_polynomial_coefficients[j];
              }

            jacobian_vector[i] = magnetic_permeability_derivative;
          }
      }
    else
      {
        std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
      }
  }

private:
  const std::vector<double> magnetic_permeability_polynomial_coefficients;
};

#endif
