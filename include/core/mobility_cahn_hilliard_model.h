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

#ifndef lethe_mobility_cahn_hilliard_model_h
#define lethe_mobility_cahn_hilliard_model_h

#include <core/interface_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * mobility for the Cahn-Hilliard-Navier-Stokes equations.
 */
class MobilityCahnHilliardModel : public InterfacePropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a MobilityCahnHilliardModel object by casting it to
   * the proper child class
   *
   * @param material_interaction_parameters Parameters for the mobility calculation
   */
  static std::shared_ptr<MobilityCahnHilliardModel>
  model_cast(
    const Parameters::MaterialInteractions &material_interaction_parameters);

  /**
   * @brief Pure virtual method to access the mobility constant
   * @return value of the mobility constant
   */
  virtual double
  get_mobility_constant() = 0;

  /**
   * @brief Definition of a virtual destructor
   */
  virtual ~MobilityCahnHilliardModel() = default;
};

/**
 * @brief Constant mobility_cahn_hilliard_constant.
 */
class MobilityCahnHilliardModelConstant : public MobilityCahnHilliardModel
{
public:
  /**
   * @brief Default constructor
   */
  MobilityCahnHilliardModelConstant(
    const double p_mobility_cahn_hilliard_constant)
    : mobility_cahn_hilliard_constant(p_mobility_cahn_hilliard_constant)
  {}

  /**
   * @brief Destructor of derived class
   */
  ~MobilityCahnHilliardModelConstant() = default;

  /**
   * @brief Method to access the mobility constant, though it returns the same value as the value function, it is implemented here for generality.
   * @return value of the mobility constant
   */
  double
  get_mobility_constant() override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief value Computes mobility_cahn_hilliard.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief vector_value Calculates the vector of mobility_cahn_hilliard.
   * @param field_vectors Vectors of the fields on which the mobility_cahn_hilliard may depend.
   * @param property_vector Vectors of the mobility_cahn_hilliard values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              mobility_cahn_hilliard_constant);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the density with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the density with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the density with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the density with respect to the field id.
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double mobility_cahn_hilliard_constant;
};

/**
 * @brief Quartic mobility_cahn_hilliard.
 */
class MobilityCahnHilliardModelQuartic : public MobilityCahnHilliardModel
{
public:
  /**
   * @brief Default constructor
   */
  MobilityCahnHilliardModelQuartic(
    const double p_mobility_cahn_hilliard_constant)
    : mobility_cahn_hilliard_constant(p_mobility_cahn_hilliard_constant)
  {
    this->model_depends_on[field::phase_order_cahn_hilliard] = true;
  }

  /**
   * @brief Destructor of derived class
   */
  ~MobilityCahnHilliardModelQuartic() = default;

  /**
   * @brief Method to access the mobility constant, though it returns the same value as the value function, it is implemented here for generality
   * @return value of the mobility constant
   */
  double
  get_mobility_constant() override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief value Calculates the mobility_cahn_hilliard.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    const double &phase_order_cahn_hilliard =
      fields_value.at(field::phase_order_cahn_hilliard);

    // The phase order values are clamped to avoid unphysical mobilities in the
    // bulk phases.
    return mobility_cahn_hilliard_constant *
           std::pow((1 - std::min(1.0,
                                  phase_order_cahn_hilliard *
                                    phase_order_cahn_hilliard)),
                    2);
  }

  /**
   * @brief vector_value Calculates the vector of mobility_cahn_hilliard_constant.
   * @param field_vectors Vectors of the fields on which the mobility_cahn_hilliard_constant may depend.
   * @param property_vector Vectors of the mobility_cahn_hilliard_constant values
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &phase_order_cahn_hilliard =
      field_vectors.at(field::phase_order_cahn_hilliard);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      {
        property_vector[i] =
          mobility_cahn_hilliard_constant *
          std::pow((1 - std::min(1.0,
                                 phase_order_cahn_hilliard[i] *
                                   phase_order_cahn_hilliard[i])),
                   2);
      }
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the density with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the density with respect to the field.
   */

  double
  jacobian(const std::map<field, double> &fields_value, field /*id*/) override
  {
    const double &phase_order_cahn_hilliard =
      fields_value.at(field::phase_order_cahn_hilliard);

    return -4 *
           (std::max(-1.0, std::min(phase_order_cahn_hilliard, 0.0)) +
            std::min(1.0, std::max(phase_order_cahn_hilliard, 0.0))) *
           mobility_cahn_hilliard_constant *
           (1 -
            std::min(1.0,
                     phase_order_cahn_hilliard * phase_order_cahn_hilliard));
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the density with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the density with respect to the field id.
   */

  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field /*id*/,
                  std::vector<double> &jacobian_vector) override
  {
    const std::vector<double> &phase_order_cahn_hilliard =
      field_vectors.at(field::phase_order_cahn_hilliard);
    for (unsigned int i = 0; i < jacobian_vector.size(); ++i)
      jacobian_vector[i] =
        -4 *
        (std::max(-1.0, std::min(phase_order_cahn_hilliard[i], 0.0)) +
         std::min(1.0, std::max(phase_order_cahn_hilliard[i], 0.0))) *
        mobility_cahn_hilliard_constant *
        (1 -
         std::min(1.0,
                  phase_order_cahn_hilliard[i] * phase_order_cahn_hilliard[i]));
  }

private:
  const double mobility_cahn_hilliard_constant;
};

#endif
