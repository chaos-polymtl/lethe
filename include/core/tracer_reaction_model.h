// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_reaction_prefactor_model_h
#define lethe_tracer_reaction_prefactor_model_h

#include <core/physical_property_model.h>

/**
 * @brief Abstract class to calculate the tracer reaction prefactor.
 *
 * This model computes the prefactor
 * \f$ k = \alpha \, C^{(n-1)} \f$
 * such that the full reaction term is assembled as
 * \f$ r = -k \, C \f$.
 *
 * Here, \f$ \alpha \f$ represents the reaction constant.
 * This parent class could have children classes that are not based on the power
 * law.
 */
class TracerReactionPrefactorModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a TracerReactionPrefactorModel object by casting it to
   * the proper child class.
   *
   * @param material_properties Parameters for a material.
   */
  static std::shared_ptr<TracerReactionPrefactorModel>
  model_cast(const Parameters::Material &material_properties);
};

/**
 * @brief Dummy tracer reaction prefactor model.
 *
 * This class is a dummy implementation and always returns zero.
 */
class NoneTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Default constructor.
   */
  NoneTracerReactionPrefactor()
  {}

  /**
   * @brief Returns a dummy value.
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return 0.;
  }

  /**
   * @brief Dummy vector value computation.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> & /*property_vector*/) override
  {}

  /**
   * @brief Returns a dummy jacobian value.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief Dummy vector jacobian computation.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> & /*jacobian_vector*/) override
  {}
};

/**
 * @brief Constant tracer reaction prefactor.
 *
 * This model computes the prefactor as
 * \f[
 *   k = \alpha \, C^{(n-1)},
 * \f]
 * where \f$ \alpha \f$ (the reaction constant) and the reaction order \f$ n \f$
 * are provided.
 */
class ConstantTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Constructor.
   *
   * @param p_tracer_reaction_constant The base reaction constant \f$ \alpha \f$.
   * @param p_tracer_reaction_order The reaction order \f$ n \f$.
   */
  ConstantTracerReactionPrefactor(const double p_tracer_reaction_constant,
                                  const double p_tracer_reaction_order)
    : tracer_reaction_constant(p_tracer_reaction_constant)
    , tracer_reaction_order(p_tracer_reaction_order)
  {
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Compute the tracer reaction prefactor.
   *
   * @param fields_value Map of field values.
   * @return The computed prefactor value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    Assert(fields_value.find(field::tracer_concentration) != fields_value.end(),
           PhysicialPropertyModelFieldUndefined(
             "ConstantTracerReactionPrefactor", "tracer_concentration"));
    return tracer_reaction_constant *
           pow(fields_value.at(field::tracer_concentration),
               tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute the vector of tracer reaction prefactor values.
   *
   * @param field_vectors Map of field vectors.
   * @param property_vector Vector to be filled with computed prefactor values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionPrefactor",
                                           "tracer_concentration"));
    const auto &concentration_vector =
      field_vectors.at(field::tracer_concentration);
    for (size_t i = 0; i < property_vector.size(); ++i)
      {
        property_vector[i] =
          tracer_reaction_constant *
          std::pow(concentration_vector[i], tracer_reaction_order - 1.);
      }
  }

  /**
   * @brief Compute the jacobian (partial derivative) of the prefactor with respect to a field.
   *
   * @param field_values Map of field values.
   * @param id The field identifier with respect to which the derivative is computed.
   * @return The computed derivative value.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "ConstantTracerReactionPrefactor", "tracer_concentration"));
    if (id == field::tracer_concentration)
      return tracer_reaction_constant * (tracer_reaction_order - 1.) *
             pow(field_values.at(field::tracer_concentration),
                 tracer_reaction_order - 2.);
    else
      return 0;
  }

  /**
   * @brief Compute the vector jacobian (partial derivatives) of the prefactor with respect to a field.
   *
   * @param field_vectors Map of field vectors.
   * @param id The field identifier with respect to which the derivative is computed.
   * @param jacobian_vector Vector to be filled with computed derivative values.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionPrefactor",
                                           "tracer_concentration"));
    if (id == field::tracer_concentration)
      {
        const auto &concentration_vector =
          field_vectors.at(field::tracer_concentration);
        for (size_t i = 0; i < jacobian_vector.size(); ++i)
          {
            jacobian_vector[i] =
              tracer_reaction_constant * (tracer_reaction_order - 1.) *
              std::pow(concentration_vector[i], tracer_reaction_order - 2.);
          }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double tracer_reaction_constant;
  const double tracer_reaction_order;
};

/**
 * @brief Tanh level-set-dependent tracer reaction prefactor.
 *
 * This model computes a prefactor that depends on the level-set function,
 * smoothly transitioning between a reaction constant \(\alpha\) inside and
 * outside the solid.
 */
class TanhLevelsetTracerReactionPrefactor : public TracerReactionPrefactorModel
{
public:
  /**
   * @brief Constructor.
   *
   * @param p_tracer_reaction_constant_outside Reaction constant \(\alpha\) outside the solid.
   * @param p_tracer_reaction_constant_inside Reaction constant \(\alpha\) inside the solid.
   * @param p_thickness Thickness for the tanh smoothing function.
   * @param p_tracer_reaction_order The reaction order.
   */
  TanhLevelsetTracerReactionPrefactor(
    const double p_tracer_reaction_constant_outside,
    const double p_tracer_reaction_constant_inside,
    const double p_thickness,
    const double p_tracer_reaction_order)
    : tracer_reaction_constant_outside(p_tracer_reaction_constant_outside)
    , tracer_reaction_constant_inside(p_tracer_reaction_constant_inside)
    , thickness(p_thickness)
    , delta_reaction_constant(tracer_reaction_constant_outside -
                              tracer_reaction_constant_inside)
    , tracer_reaction_order(p_tracer_reaction_order)
  {
    this->model_depends_on[field::levelset]             = true;
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Compute the tracer reaction prefactor.
   *
   * @param field_values Map of field values.
   * @return The computed prefactor value.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    const double levelset      = field_values.at(field::levelset);
    const double concentration = field_values.at(field::tracer_concentration);
    const double k =
      tracer_reaction_constant_inside +
      delta_reaction_constant * (0.5 + 0.5 * tanh(levelset / thickness));
    return k * pow(concentration, tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute a vector of tracer reaction prefactor values.
   *
   * @param field_vectors Map of field vectors.
   * @param property_vector Vector to be filled with computed prefactor values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));

    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);
    const std::vector<double> &concentration_vec =
      field_vectors.at(field::tracer_concentration);

    const unsigned int n_values = levelset_vec.size();
    Assert(n_values == levelset_vec.size(),
           SizeOfFields(n_values, levelset_vec.size()));
    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double k = tracer_reaction_constant_inside +
                         delta_reaction_constant *
                           (0.5 + 0.5 * tanh(levelset_vec[i] / thickness));
        property_vector[i] =
          k * pow(concentration_vec[i], tracer_reaction_order - 1.);
      }
  }

  /**
   * @brief Compute the jacobian (partial derivative) of the prefactor with respect to a specified field.
   *
   * @param field_values Map of field values.
   * @param id The field identifier with respect to which the derivative is computed.
   * @return The computed derivative value.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_values.find(field::tracer_concentration) != field_values.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    if (id == field::levelset)
      {
        return numerical_jacobian(field_values, field::levelset);
      }
    else if (id == field::tracer_concentration)
      {
        const double k =
          tracer_reaction_constant_inside +
          delta_reaction_constant *
            (0.5 + 0.5 * tanh(field_values.at(field::levelset) / thickness));
        return k * (tracer_reaction_order - 1.) *
               pow(field_values.at(field::tracer_concentration),
                   tracer_reaction_order - 2.);
      }
    else
      return 0;
  }

  /**
   * @brief Compute the vector jacobian (partial derivatives) of the prefactor with respect to a field.
   *
   * @param field_vectors Map of field vectors.
   * @param id The field identifier with respect to which the derivative is computed.
   * @param jacobian_vector Vector to be filled with computed derivative values.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "levelset"));
    Assert(field_vectors.find(field::tracer_concentration) !=
             field_vectors.end(),
           PhysicialPropertyModelFieldUndefined(
             "TanhLevelsetTracerReactionPrefactor", "tracer_concentration"));
    if (id == field::levelset)
      {
        vector_numerical_jacobian(field_vectors, id, jacobian_vector);
      }
    else if (id == field::tracer_concentration)
      {
        const std::vector<double> &levelset_vec =
          field_vectors.at(field::levelset);
        const std::vector<double> &concentration_vec =
          field_vectors.at(field::tracer_concentration);
        const unsigned int n_values = levelset_vec.size();
        for (unsigned int i = 0; i < n_values; ++i)
          {
            const double k = tracer_reaction_constant_inside +
                             delta_reaction_constant *
                               (0.5 + 0.5 * tanh(levelset_vec[i] / thickness));
            jacobian_vector[i] =
              k * (tracer_reaction_order - 1.) *
              pow(concentration_vec[i], tracer_reaction_order - 2.);
          }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.);
  }

private:
  const double tracer_reaction_constant_outside;
  const double tracer_reaction_constant_inside;
  const double thickness;
  const double delta_reaction_constant;
  const double tracer_reaction_order;
};

#endif
