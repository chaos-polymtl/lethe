// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_reaction_constant_model_h
#define lethe_tracer_reaction_constant_model_h

#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculates the
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
 * @brief Constant tracer reaction constant.
 */
class ConstantTracerReactionConstant : public TracerReactionConstantModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantTracerReactionConstant(const double p_tracer_reaction_constant)
    : tracer_reaction_constant(p_tracer_reaction_constant)
  {}

  /**
   * @brief value Calculates the tracer reaction constant
   * @param fields_value Value of the various fields on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return tracer_reaction_constant;
  };

  /**
   * @brief vector_value Calculates the vector of tracer reaction constant.
   * @param field_vectors Vectors of the fields on which the reaction constant may depend.
   * @param property_vector Vectors of the tracer reaction constant values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              tracer_reaction_constant);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the reaction constant with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the reaction constant with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculates the derivative of the tracer reaction constant with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the tracer reaction constant with respect to the field id.
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
  const double tracer_reaction_constant;
};

/**
 * @brief Reaction constant that depends on the levelset
 */
class TanhLevelsetTracerReactionConstant : public TracerReactionConstantModel
{
public:
  /**
   * @brief Constructor of the levelset-dependent reaction constant model.
   *
   * @param[in] p_tracer_reaction_constant_outside Reaction constant outside the
   * solid
   * @param[in] p_tracer_reaction_constant_inside Reaction constant inside the
   * solid
   * @param[in] p_thickness Thickness of the tanh function used to smooth the
   * property jump
   */
  TanhLevelsetTracerReactionConstant(
    const double p_tracer_reaction_constant_outside,
    const double p_tracer_reaction_constant_inside,
    const double p_thickness)
    : tracer_reaction_constant_outside(p_tracer_reaction_constant_outside)
    , tracer_reaction_constant_inside(p_tracer_reaction_constant_inside)
    , thickness(p_thickness)
    , delta_reaction_constant(tracer_reaction_constant_outside -
                              tracer_reaction_constant_inside)
  {
    this->model_depends_on[field::levelset] = true;
  }

  /**
   * @brief Compute the reaction constant.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the reaction constant depends on the levelset.
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
    double levelset = field_values.at(field::levelset);

    return tracer_reaction_constant_inside +
           delta_reaction_constant * (0.5 + 0.5 * tanh(levelset / thickness));
  }

  /**
   * @brief Compute a vector of reaction constant.
   *
   * @param[in] field_vectors Vectors of the fields on which the reaction
   * constant may depend. In this case, the reaction constant depends on the
   * levelset. The map stores a vector of values per field.
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

    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);

    const unsigned int n_values = levelset_vec.size();

    Assert(n_values == levelset_vec.size(),
           SizeOfFields(n_values, levelset_vec.size()));

    for (unsigned int i = 0; i < n_values; ++i)
      {
        double levelset = levelset_vec[i];
        property_vector[i] =
          tracer_reaction_constant_inside +
          delta_reaction_constant * (0.5 + 0.5 * tanh(levelset / thickness));
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
    if (id == field::levelset)
      {
        AssertThrow(field_values.find(field::levelset) != field_values.end(),
                    PhysicialPropertyModelFieldUndefined(
                      "TanhLevelsetTracerReactionConstant", "levelset"));
        return numerical_jacobian(field_values, field::levelset);
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
    if (id == field::levelset)
      {
        AssertThrow(field_vectors.find(field::levelset) != field_vectors.end(),
                    PhysicialPropertyModelFieldUndefined(
                      "TanhLevelsetTracerReactionConstant", "levelset"));
        vector_numerical_jacobian(field_vectors, id, jacobian_vector);
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.);
  };


private:
  const double tracer_reaction_constant_outside;
  const double tracer_reaction_constant_inside;
  const double thickness;
  const double delta_reaction_constant;
};

#endif
