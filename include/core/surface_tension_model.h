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

#ifndef lethe_surface_tension_model_h
#define lethe_surface_tension_model_h

#include <core/interface_property_model.h>

/**
 * @brief SurfaceTensionModel. Abstract class that allows to calculate the
 * surface tension coefficient.
 */
class SurfaceTensionModel : public InterfacePropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a SurfaceTensionModel object
   * by casting it to the proper child class
   *
   * @param material_interaction_parameters Parameters for the surface tension
   * coefficient calculation
   */
  static std::shared_ptr<SurfaceTensionModel>
  model_cast(
    const Parameters::MaterialInteractions &material_interaction_parameters);

  /**
   * @brief is_constant_surface_tension_model Returns a boolean indicating if
   * the model is a constant surface tension model.
   * @return Boolean value of if the model corresponds to a constant surface
   * tension model.
   */
  bool
  is_constant_surface_tension_model()
  {
    return surface_tension_is_constant;
  }

protected:
  bool surface_tension_is_constant = false;
};


/**
 * @brief Constant surface tension.
 */
class SurfaceTensionConstant : public SurfaceTensionModel
{
public:
  /**
   * @brief Default constructor
   */
  SurfaceTensionConstant(const double p_surface_tension_coefficient)
    : surface_tension_coefficient(p_surface_tension_coefficient)
  {
    this->surface_tension_is_constant = true;
  }

  /**
   * @brief value Calculates the surface tension coefficient
   * @param fields_value Value of the various field on which the property may
   * depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return surface_tension_coefficient;
  }

  /**
   * @brief vector_value Calculates the vector of surface tension coefficient.
   * @param field_vectors Vectors of the fields on which the surface tension
   * coefficient may depend.
   * @param property_vector Vectors of the surface tension coefficient values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              surface_tension_coefficient);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * surface tension coefficient with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the surface tension coefficient
   * with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the surface tension
   * coefficient with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian vector of the value of the derivative of the surface
   * tension coefficient with respect to the field id.
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
  const double surface_tension_coefficient;
};

/**
 * @brief Linear surface tension. The surface tension is given by:
 * sigma_0 + dsimga/dT * (T-T_0).
 */
class SurfaceTensionLinear : public SurfaceTensionModel
{
public:
  /**
   * @brief Default constructor
   */
  SurfaceTensionLinear(
    const Parameters::SurfaceTensionParameters &p_surface_tension_parameters)
    : surface_tension_coefficient(
        p_surface_tension_parameters.surface_tension_coefficient)
    , T_0(p_surface_tension_parameters.T_0)
    , surface_tension_gradient(
        p_surface_tension_parameters.surface_tension_gradient)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief value Calculates the surface tension coefficient
   * @param field_vectors Vectors of the fields on which the surface tension
   * depends. In this case, it depends on the temperature.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    const double temperature = fields_value.at(field::temperature);
    return surface_tension_coefficient +
           surface_tension_gradient * (temperature - T_0);
  }

  /**
   * @brief vector_value Calculates the vector of surface tension coefficient.
   * @param field_vectors Vectors of the fields on which the surface tension
   * depends. In this case, it depends on the temperature.
   * @param property_vector Vectors of the surface tension coefficient values
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      property_vector[i] = surface_tension_coefficient +
                           surface_tension_gradient * (temperature[i] - T_0);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * surface tension coefficient with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the surface tension coefficient
   * with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/, field id) override
  {
    if (id == field::temperature)
      return surface_tension_gradient;
    else
      return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the surface tension
   * coefficient with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property.
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated.
   * @param jacobian vector of the value of the derivative of the surface
   * tension coefficient with respect to the field id.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field          id,
    std::vector<double> &jacobian_vector) override
  {
    if (id == field::temperature)
      std::fill(jacobian_vector.begin(),
                jacobian_vector.end(),
                surface_tension_gradient);
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double surface_tension_coefficient;
  const double T_0;
  const double surface_tension_gradient;
};

#endif
