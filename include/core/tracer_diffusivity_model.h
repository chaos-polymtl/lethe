// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_diffusivity_model_h
#define lethe_tracer_diffusivity_model_h

#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculates the
 * tracer diffusivity.
 */
class TracerDiffusivityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a TracerDiffusivityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<TracerDiffusivityModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant tracer diffusivity.
 */
class ConstantTracerDiffusivity : public TracerDiffusivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantTracerDiffusivity(const double p_tracer_diffusivity)
    : tracer_diffusivity(p_tracer_diffusivity)
  {}

  /**
   * @brief value Calculates the tracer diffusivity
   * @param fields_value Value of the various fields on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return tracer_diffusivity;
  };

  /**
   * @brief vector_value Calculates the vector of tracer diffusivity.
   * @param field_vectors Vectors of the fields on which the diffusivity may depend.
   * @param property_vector Vectors of the tracer diffusivity values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              tracer_diffusivity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the diffusivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the diffusivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculates the derivative of the tracer diffusivity with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the tracer diffusivity with respect to the field id.
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
  const double tracer_diffusivity;
};

/**
 * @brief Diffusivity that depends on the level set
 */
class TanhLevelsetTracerDiffusivity : public TracerDiffusivityModel
{
public:
  /**
   * @brief Constructor of the level set-dependent diffusivity model.
   *
   * @param[in] p_tracer_diffusivity_outside Diffusivity outside the solid
   * @param[in] p_tracer_diffusivity_inside Diffusivity inside the solid
   * @param[in] p_thickness Thickness of the tanh function used to smooth the
   * property jump
   */
  TanhLevelsetTracerDiffusivity(const double p_tracer_diffusivity_outside,
                                const double p_tracer_diffusivity_inside,
                                const double p_thickness)
    : tracer_diffusivity_outside(p_tracer_diffusivity_outside)
    , tracer_diffusivity_inside(p_tracer_diffusivity_inside)
    , thickness(p_thickness)
    , delta_diffusivity(tracer_diffusivity_outside - tracer_diffusivity_inside)
  {
    this->model_depends_on[field::levelset] = true;
  }

  /**
   * @brief Compute the diffusivity.
   *
   * @param[in] field_values Values of the various fields on which the property
   * may depend. In this case, the diffusivity depends on the level set.
   * The map stores a single value per field.
   *
   * @return Value of the diffusivity computed with the @p field_values.
   */
  double
  value(const std::map<field, double> &field_values) override
  {
    Assert(field_values.find(field::levelset) != field_values.end(),
           PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerDiffusivity",
                                                "levelset"));
    double levelset = field_values.at(field::levelset);

    return tracer_diffusivity_inside +
           delta_diffusivity * (0.5 + 0.5 * tanh(levelset / thickness));
  }

  /**
   * @brief Compute a vector of diffusivity.
   *
   * @param[in] field_vectors Vectors of the fields on which the diffusivity
   * may depend. In this case, the diffusivity depends on the level set. The map
   * stores a vector of values per field.
   *
   * @param[out] property_vector Vectors of computed diffusivities.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::levelset) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("TanhLevelsetTracerDiffusivity",
                                                "levelset"));

    const std::vector<double> &levelset_vec = field_vectors.at(field::levelset);

    const unsigned int n_values = property_vector.size();

    Assert(n_values == levelset_vec.size(),
           SizeOfFields(n_values, levelset_vec.size()));

    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double levelset = levelset_vec[i];
        property_vector[i] =
          tracer_diffusivity_inside +
          delta_diffusivity * (0.5 + 0.5 * tanh(levelset / thickness));
      }
  }

  /**
   * @brief Compute the jacobian (the partial derivative) of the diffusivity
   * with respect to a specified field.
   *
   * @param[in] field_values Values of the various fields on which the specific
   * heat may depend. The map stores a single value per field.
   *
   * @param[in] id Indicator of the field with respect to which the jacobian
   * should be computed.
   *
   * @return Value of the derivative of the diffusivity with respect to the
   * specified field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::levelset)
      {
        Assert(field_values.find(field::levelset) != field_values.end(),
               PhysicialPropertyModelFieldUndefined(
                 "TanhLevelsetTracerDiffusivity", "levelset"));
        return numerical_jacobian(field_values, field::levelset);
      }
    else
      return 0;
  };

  /**
   * @brief Compute the derivative of the diffusivity with respect to a field.
   *
   * @param[in] field_vectors Vector of values of the fields used to evaluate
   * the diffusivity. The map stores a vector of values per field.
   *
   * @param[in] id Identifier of the field with respect to which a derivative
   * should be computed.
   *
   * @param[out] jacobian Vector of computed derivative values of the
   * diffusivity with respect to the field of the specified @p id.
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    if (id == field::levelset)
      {
        Assert(field_vectors.find(field::levelset) != field_vectors.end(),
               PhysicialPropertyModelFieldUndefined(
                 "TanhLevelsetTracerDiffusivity", "levelset"));
        vector_numerical_jacobian(field_vectors, id, jacobian_vector);
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0.);
  };


private:
  const double tracer_diffusivity_outside;
  const double tracer_diffusivity_inside;
  const double thickness;
  const double delta_diffusivity;
};

#endif
