/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_thermal_conductivity_model_h
#define lethe_thermal_conductivity_model_h

#include <core/parameters.h>
#include <core/phase_change.h>
#include <core/physical_property_model.h>

/**
 * @brief ThermalConductivityModel. Abstract class that allows to calculate the
 * thermal conductivity on each quadrature point using the temperature of the
 * fluid.
 */
class ThermalConductivityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ThermalConductivityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ThermalConductivityModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant thermal conductivity.
 */
class ConstantThermalConductivity : public ThermalConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantThermalConductivity(const double p_thermal_conductivity)
    : thermal_conductivity(p_thermal_conductivity)
  {}

  /**
   * @brief value Calculates the value the thermal conductivity
   * @param fields_value Value of the various field on which the thermal conductivity depends.
   * @return value of the thermal conductivity calculated with the fields_value.
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return thermal_conductivity;
  };

  /**
   * @brief vector_value Calculates the vector value of thermal conductivities
   * @param field_vectors Vector of properties on which the thermal conductivities depend
   * @param property_vector Values of the thermal conductivities
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    property_vector.assign(property_vector.size(), thermal_conductivity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the thermal conductivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the thermal conductivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the thermal conductivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the thermal conductivity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double thermal_conductivity;
};


/**
 * @brief ThermalConductivityLinear Implements a linear temperature-dependant thermal conductivity of the form k = A + BT
 */
class ThermalConductivityLinear : public ThermalConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ThermalConductivityLinear(const double A, const double B)
    : A(A)
    , B(B)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief value Calculates the value the thermal conductivity
   * @param fields_value Value of the various field on which the thermal conductivity depends.
   * @return value of the thermal conductivity calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    return A + B * fields_value.at(field::temperature);
  };

  /**
   * @brief vector_value Calculates the vector value of thermal conductivities
   * @param field_vectors Vector of properties on which the thermal conductivities depend
   * @param property_vector Values of the thermal conductivities
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &T = field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      property_vector[i] = A + B * T[i];
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the thermal conductivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the thermal conductivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return B;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the thermal conductivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the thermal conductivity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), B);
  };

private:
  const double A;
  const double B;
};

/**
 * @brief ThermalConductivityPhaseChange Implements a phase-dependant thermal conductivity
 */
class ThermalConductivityPhaseChange : public ThermalConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ThermalConductivityPhaseChange(
    const Parameters::PhaseChange phase_change_params)
    : p_phase_change_params(phase_change_params)
    , thermal_conductivity_s(phase_change_params.thermal_conductivity_s)
    , thermal_conductivity_l(phase_change_params.thermal_conductivity_l)
    , T_solidus(phase_change_params.T_solidus)
    , T_liquidus(phase_change_params.T_liquidus)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief value Calculates the value of the thermal conductivity
   * @param fields_value Value of the various field on which the thermal conductivity depends.
   * @return value of the thermal conductivity calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    double thermal_conductivity;
    // Thermal conductivity of solid phase
    if (fields_value.at(field::temperature) < T_solidus)
      thermal_conductivity = thermal_conductivity_s;
    // Thermal conductivity of liquid phase
    else if (fields_value.at(field::temperature) > T_liquidus)
      thermal_conductivity = thermal_conductivity_l;
    // Mean value of the conductivities of the solid and liquid phases
    else
      {
        const double l_frac =
          calculate_liquid_fraction(fields_value.at(field::temperature),
                                    p_phase_change_params);

        thermal_conductivity = thermal_conductivity_l * l_frac +
                               thermal_conductivity_s * (1. - l_frac);
      }

    return thermal_conductivity;
  };

  /**
   * @brief vector_value Calculates the vector value of thermal conductivities
   * @param field_vectors Vector of properties on which the thermal conductivities depend
   * @param property_vector Values of the thermal conductivities
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &T = field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      {
        // Thermal conductivity of solid phase
        if (T[i] < T_solidus)
          property_vector[i] = thermal_conductivity_s;
        // Thermal conductivity of liquid phase
        else if (T[i] > T_liquidus)
          property_vector[i] = thermal_conductivity_l;
        else
          {
            const double l_frac =
              calculate_liquid_fraction(T[i], p_phase_change_params);

            property_vector[i] = thermal_conductivity_l * l_frac +
                                 thermal_conductivity_s * (1. - l_frac);
          }
      }
  };

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the thermal conductivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the thermal conductivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::temperature)
      return this->numerical_jacobian(field_values, field::temperature);
    else
      return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the thermal conductivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the thermal conductivity with respect to the field id
   */

  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    vector_numerical_jacobian(field_vectors, id, jacobian_vector);
  };

private:
  Parameters::PhaseChange p_phase_change_params;
  const double            thermal_conductivity_s;
  const double            thermal_conductivity_l;
  const double            T_solidus;
  const double            T_liquidus;
};

#endif
