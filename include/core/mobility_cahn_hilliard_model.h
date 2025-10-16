// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_mobility_cahn_hilliard_model_h
#define lethe_mobility_cahn_hilliard_model_h

#include <core/interface_property_model.h>

/**
 * @brief Implementation of the computation of the mobility
 * for Cahn-Hilliard equations.
 */
class MobilityCahnHilliardModel : public InterfacePropertyModel
{
public:
  /**
   * @brief Instantiate and return a pointer to a MobilityCahnHilliardModel
   * object by casting it to the proper child class.
   *
   * @param[in] material_interaction_parameters Parameters for the mobility
   * calculation.
   */
  static std::shared_ptr<MobilityCahnHilliardModel>
  model_cast(
    const Parameters::MaterialInteractions &material_interaction_parameters);

  /**
   * @brief Pure virtual method to access the mobility constant.
   * @return Value of the mobility constant.
   */
  virtual double
  get_mobility_constant() = 0;

  /**
   * @brief Definition of a virtual destructor for the class.
   */
  ~MobilityCahnHilliardModel() override = default;
};

/**
 * @brief Constant mobility model.
 *
 * The mobility function is the following : \f$M(\phi) = D \f$ where D is the
 * mobility constant.
 */
class MobilityCahnHilliardModelConstant : public MobilityCahnHilliardModel
{
public:
  /**
   * @brief Default constructor
   *
   * @param[in] p_mobility_cahn_hilliard_constant the user defined mobility
   * constant.
   */
  MobilityCahnHilliardModelConstant(
    const double p_mobility_cahn_hilliard_constant)
    : mobility_cahn_hilliard_constant(p_mobility_cahn_hilliard_constant)
  {}

  /**
   * @brief Destructor of derived class.
   */
  ~MobilityCahnHilliardModelConstant() override = default;

  /**
   * @brief Return the mobility constant, though it returns the same
   * value as the value method, it is implemented here for generality.
   *
   * @return Value of the mobility constant.
   */
  double
  get_mobility_constant() override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief Compute the mobility.
   * @param[in] fields_value Value of the various field on which the mobility
   * may depend.
   * @return Value of the mobility.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief Calculate the vector of mobility_cahn_hilliard.
   * @param[in] field_vectors Vectors of the fields on which the mobility
   * may depend.
   * @param[out] property_vector Vectors of the mobility values
   */
  void
  vector_value(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double> &property_vector) override
  {
    std::ranges::fill(property_vector, mobility_cahn_hilliard_constant);
  }

  /**
   * @brief Calculate the jacobian (the partial derivative) of the
   * mobility with respect to a field.
   * @param[in] field_values Value of the various fields on which the mobility
   * may depend.
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return Value of the partial derivative of the mobility with respect to
   * the field.
   */
  double
  jacobian([[maybe_unused]] const std::map<field, double> &field_values,
           [[maybe_unused]] field                          id) override
  {
    return 0;
  }

  /**
   * @brief Calculate the derivative of the mobility with respect to a field.
   * @param[in] field_vectors Vector for the values of the fields used to
   * evaluate the mobility.
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be calculated.
   * @param[out] jacobian_vector Vector of the value of the derivative of the
   * mobility with respect to the field id.
   */
  void
  vector_jacobian(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    [[maybe_unused]] const field                                 id,
    std::vector<double> &jacobian_vector) override
  {
    std::ranges::fill(jacobian_vector, 0);
  }

private:
  const double mobility_cahn_hilliard_constant;
};

/**
 * @brief Quartic mobility model.
 *
 * The mobility function is the following : \f$M(\phi) = D(1-\phi^2)^2 \f$
 * where D is the mobility constant.
 */
class MobilityCahnHilliardModelQuartic : public MobilityCahnHilliardModel
{
public:
  /**
   * @brief Default constructor.
   *
   * @param[in] p_mobility_cahn_hilliard_constant the user defined mobility
   * constant.
   */
  MobilityCahnHilliardModelQuartic(
    const double p_mobility_cahn_hilliard_constant)
    : mobility_cahn_hilliard_constant(p_mobility_cahn_hilliard_constant)
  {
    this->model_depends_on[field::phase_order_cahn_hilliard] = true;
  }

  /**
   * @brief Destructor of derived class.
   */
  ~MobilityCahnHilliardModelQuartic() override = default;

  /**
   * @brief Return the mobility constant.
   * @return Value of the mobility constant.
   */
  double
  get_mobility_constant() override
  {
    return mobility_cahn_hilliard_constant;
  }

  /**
   * @brief Calculate the mobility.
   * @param[in] fields_value Value of the various fields on which the mobility
   * may depend.
   * @return Value of the mobility calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    Assert(fields_value.contains(field::phase_order_cahn_hilliard),
           PhysicialPropertyModelFieldUndefined(
             "MobilityCahnHilliardModelQuartic", "phase_order_cahn_hilliard"));
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
   * @brief Calculate the vector of mobility.
   * @param[in] field_vectors Vectors of the fields on which the mobility
   * may depend.
   * @param[out] property_vector Vector of the mobility values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.contains(field::phase_order_cahn_hilliard),
           PhysicialPropertyModelFieldUndefined(
             "MobilityCahnHilliardModelQuartic", "phase_order_cahn_hilliard"));
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
   * @brief Calculate the jacobian (the partial derivative) of the
   * mobility with respect to a field.
   * @param[in] fields_value Value of the various fields on which the mobility
   * may depend.
   * @param[in] id Indicator of the field with respect to which the Jacobian
   * should be calculated.
   * @return Value of the partial derivative of the mobility with respect to
   * the field.
   */

  double
  jacobian(const std::map<field, double> &fields_value,
           [[maybe_unused]] field         id) override
  {
    Assert(fields_value.contains(field::phase_order_cahn_hilliard),
           PhysicialPropertyModelFieldUndefined(
             "MobilityCahnHilliardModelQuartic", "phase_order_cahn_hilliard"));
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
   * @brief Calculate the derivative of the mobility with respect to a field.
   * @param[in] fields_vectors Vector for the values of the fields used to
   * evaluate the mobility.
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be calculated.
   * @param[out] jacobian_vector Vector of the values of the derivatives of the
   * mobility with respect to the field id.
   */

  void
  vector_jacobian(const std::map<field, std::vector<double>> &fields_vectors,
                  [[maybe_unused]] const field                id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(fields_vectors.contains(field::phase_order_cahn_hilliard),
           PhysicialPropertyModelFieldUndefined(
             "MobilityCahnHilliardModelQuartic", "phase_order_cahn_hilliard"));
    const std::vector<double> &phase_order_cahn_hilliard =
      fields_vectors.at(field::phase_order_cahn_hilliard);
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
