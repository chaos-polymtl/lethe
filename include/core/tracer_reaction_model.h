// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_reaction_constant_model_h
#define lethe_tracer_reaction_constant_model_h

#include <core/physical_property_model.h>

/**
 * @brief Abstract class to calculate the
 * tracer reaction constant.
 */
class TracerReactionConstantModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a TracerReactionConstantModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<TracerReactionConstantModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief None tracer reaction constant mode. This class is dummy.
 */
class NoneTracerReactionConstant : public TracerReactionConstantModel
{
public:
  /**
   * @brief Default constructor
   */
  NoneTracerReactionConstant()
  {}

  /**
   * @brief Function of dummy class.
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return 0.;
  };

  /**
   * @brief Function of dummy class.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> & /*property_vector*/) override
  {}

  /**
   * @brief Function of dummy class.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief Function of dummy class.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> & /*jacobian_vector*/) override{};
};

/**
 * @brief Constant tracer reaction constant.
 * // TODO Change the comments and documentation
 */
class ConstantTracerReactionConstant : public TracerReactionConstantModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantTracerReactionConstant(const double p_tracer_reaction_constant,
                                 const double p_tracer_reaction_order)
    : tracer_reaction_constant(p_tracer_reaction_constant)
    , tracer_reaction_order(p_tracer_reaction_order)
  {
    this->model_depends_on[field::tracer_concentration] = true;
  }

  /**
   * @brief Calculate the tracer reaction constant
   * @param fields_value Value of the various fields on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    // TODO change the calculation
    AssertThrow(
      fields_value.find(field::tracer_concentration) != fields_value.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionConstant",
                                           "tracer_concentration"));
    return tracer_reaction_constant *
           pow(fields_value.at(field::tracer_concentration),
               tracer_reaction_order - 1.);
  };

  /**
   * @brief Calculate the vector of tracer reaction constant.
   * @param field_vectors Vectors of the fields on which the reaction constant may depend.
   * @param property_vector Vectors of the tracer reaction constant values
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    AssertThrow(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionConstant",
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
   * @brief Calculate the jacobian (the partial derivative) of the reaction constant with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the reaction constant with respect to the field.
   */

  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    AssertThrow(
      field_values.find(field::tracer_concentration) != field_values.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionConstant",
                                           "tracer_concentration"));
    if (id == field::tracer_concentration)
      return tracer_reaction_constant *
             pow(field_values.at(field::tracer_concentration),
                 tracer_reaction_order - 2.);
    else
      return 0;
  };

  /**
   * @brief Calculate the derivative of the tracer reaction constant with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the tracer reaction constant with respect to the field id.
   */

  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    AssertThrow(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("ConstantTracerReactionConstant",
                                           "tracer_concentration"));
    if (id == field::tracer_concentration)
      {
        const auto &concentration_vector =
          field_vectors.at(field::tracer_concentration);
        for (size_t i = 0; i < jacobian_vector.size(); ++i)
          {
            jacobian_vector[i] =
              tracer_reaction_constant *
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
 * @brief Reaction constant that depends on the level set
 */
class TanhLevelsetTracerReactionConstant : public TracerReactionConstantModel
{
public:
  /**
   * @brief Constructor of the level set-dependent reaction constant model.
   *
   * @param[in] p_tracer_reaction_constant_outside Reaction constant outside the
   * solid
   * @param[in] p_tracer_reaction_constant_inside Reaction constant inside the
   * solid
   * @param[in] p_thickness Thickness of the tanh function used to smooth the
   * reaction constant jump between the inside and outside of the immersed
   * solid.
   */
  TanhLevelsetTracerReactionConstant(
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
   * @brief Compute the reaction constant.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the reaction constant depends on the level set.
   * The map stores a single value per field.
   *
   * @return Value of the reaction constant computed with the @p field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    AssertThrow(field_values.find(field::levelset) != field_values.end(),
                PhysicialPropertyModelFieldUndefined(
                  "TanhLevelsetTracerReactionConstant", "levelset"));
    AssertThrow(
      field_values.find(field::tracer_concentration) != field_values.end(),
      PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerReactionConstant",
                                           "tracer_concentration"));
    const double levelset      = field_values.at(field::levelset);
    const double concentration = field_values.at(field::tracer_concentration);
    const double k =
      tracer_reaction_constant_inside +
      delta_reaction_constant * (0.5 + 0.5 * tanh(levelset / thickness));
    return k * pow(concentration, tracer_reaction_order - 1.);
  }

  /**
   * @brief Compute a vector of reaction constant.
   *
   * @param[in] field_vectors Vectors of the fields on which the reaction
   * constant may depend. In this case, the reaction constant depends on the
   * level set. The map stores a vector of values per field.
   *
   * @param[out] property_vector Vectors of computed reaction constants.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    AssertThrow(field_vectors.find(field::levelset) != field_vectors.end(),
                PhysicialPropertyModelFieldUndefined(
                  "TanhLevelsetTracerReactionConstant", "levelset"));
    AssertThrow(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerReactionConstant",
                                           "tracer_concentration"));

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
   * @brief Compute the jacobian (the partial derivative) of the reaction constant
   * with respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @return Value of the derivative of the reaction constant with respect to the
   * specified field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    AssertThrow(field_values.find(field::levelset) != field_values.end(),
                PhysicialPropertyModelFieldUndefined(
                  "TanhLevelsetTracerReactionConstant", "levelset"));
    AssertThrow(
      field_values.find(field::tracer_concentration) != field_values.end(),
      PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerReactionConstant",
                                           "tracer_concentration"));
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
        return k * pow(field_values.at(field::tracer_concentration),
                       tracer_reaction_order - 2.);
      }
    else
      return 0;
  };

  /**
   * @brief Compute the derivative of the reaction constant with respect to a field.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the reaction constant. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian Vector of computed derivative values of the
   * reaction constant with respect to the field of the specified @p id.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    AssertThrow(field_vectors.find(field::levelset) != field_vectors.end(),
                PhysicialPropertyModelFieldUndefined(
                  "TanhLevelsetTracerReactionConstant", "levelset"));
    AssertThrow(
      field_vectors.find(field::tracer_concentration) != field_vectors.end(),
      PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerReactionConstant",
                                           "tracer_concentration"));
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
              k * pow(concentration_vec[i], tracer_reaction_order - 2.);
          }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.);
  };


private:
  const double tracer_reaction_constant_outside;
  const double tracer_reaction_constant_inside;
  const double thickness;
  const double delta_reaction_constant;
  const double tracer_reaction_order;
};

#endif
