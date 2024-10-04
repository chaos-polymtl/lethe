// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_surface_tension_model_h
#define lethe_surface_tension_model_h

#include <core/interface_property_model.h>
#include <core/phase_change.h>

/**
 * @brief Abstract class that allows to calculate the
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
   * @param p_surface_tension_parameters The parameters for the surface tension
   * model, including the surface tension coefficient sigma_0 and surface
   * tension gradient dsimga/dT at the reference temperature T_0
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
    AssertThrow(fields_value.find(field::temperature) != fields_value.end(),
                PhysicialPropertyModelFieldUndefined("SurfaceTensionLinear",
                                                     "temperature"));
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
    AssertThrow(field_vectors.find(field::temperature) != field_vectors.end(),
                PhysicialPropertyModelFieldUndefined("SurfaceTensionLinear",
                                                     "temperature"));
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

protected:
  const double surface_tension_coefficient;
  const double T_0;
  const double surface_tension_gradient;
};

/**
 * @brief Linear surface tension with phase change. The surface tension is given by:
 * alpha_l*(sigma_0 + dsigma/dT * (T-T_0)).
 */
class SurfaceTensionPhaseChange : public SurfaceTensionLinear
{
public:
  /**
   * @brief Default constructor
   * @param p_surface_tension_parameters The parameters for the surface tension
   * model, including the surface tension coefficient sigma_0 and surface
   * tension gradient dsimga/dT at the reference temperature T_0, and the
   * solidus and liquidus temperatures T_solidus and T_liquidus respectively.
   */
  SurfaceTensionPhaseChange(
    const Parameters::SurfaceTensionParameters &p_surface_tension_parameters)
    : SurfaceTensionLinear(p_surface_tension_parameters)
    , T_solidus(p_surface_tension_parameters.T_solidus)
    , T_liquidus(p_surface_tension_parameters.T_liquidus)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief value Calculates the surface tension coefficient.
   * @param field_vectors Vectors of the fields on which the surface tension
   * depends. In this case, it depends on the temperature.
   * @return value of the physical property calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    AssertThrow(fields_value.find(field::temperature) != fields_value.end(),
                PhysicialPropertyModelFieldUndefined(
                  "SurfaceTensionPhaseChange", "temperature"));
    const double temperature = fields_value.at(field::temperature);

    double surface_tension;

    if (temperature < T_solidus)
      surface_tension = 0.0;
    else if (temperature > T_liquidus)
      surface_tension = surface_tension_coefficient +
                        surface_tension_gradient * (temperature - T_0);
    else
      {
        const double l_frac =
          calculate_liquid_fraction(temperature, T_solidus, T_liquidus);

        surface_tension = (surface_tension_coefficient +
                           surface_tension_gradient * (temperature - T_0)) *
                          l_frac;
      }
    return surface_tension;
  }

  /**
   * @brief vector_value Calculates the vector of surface tension coefficient.
   * @param field_vectors Vectors of the fields on which the surface tension
   * depends. In this case, it depends on the temperature.
   * @param property_vector Vectors of the surface tension coefficient values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    AssertThrow(field_vectors.find(field::temperature) != field_vectors.end(),
                PhysicialPropertyModelFieldUndefined(
                  "SurfaceTensionPhaseChange", "temperature"));
    const std::vector<double> &temperature =
      field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      if (temperature[i] < T_solidus)
        property_vector[i] = 0.0;
      else if (temperature[i] > T_liquidus)
        property_vector[i] = surface_tension_coefficient +
                             surface_tension_gradient * (temperature[i] - T_0);
      else
        {
          const double l_frac =
            calculate_liquid_fraction(temperature[i], T_solidus, T_liquidus);

          property_vector[i] =
            (surface_tension_coefficient +
             surface_tension_gradient * (temperature[i] - T_0)) *
            l_frac;
        }
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * surface tension coefficient with respect to a field.
   * @param field_values Value of the various fields on which the property may
   * depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the surface tension coefficient
   * with respect to the field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::temperature)
      {
        AssertThrow(field_values.find(field::temperature) != field_values.end(),
                    PhysicialPropertyModelFieldUndefined(
                      "SurfaceTensionPhaseChange", "temperature"));
        const double temperature = field_values.at(field::temperature);
        if (temperature < T_solidus)
          return 0;
        else if (temperature > T_liquidus)
          return surface_tension_gradient;
        else
          {
            const double l_frac =
              calculate_liquid_fraction(temperature, T_solidus, T_liquidus);
            return surface_tension_gradient * l_frac;
          }
      }
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
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    if (id == field::temperature)
      {
        AssertThrow(
          field_vectors.find(field::temperature) != field_vectors.end(),
          PhysicialPropertyModelFieldUndefined("SurfaceTensionPhaseChange",
                                               "temperature"));
        const std::vector<double> &temperature =
          field_vectors.at(field::temperature);
        for (unsigned int i = 0; i < jacobian_vector.size(); ++i)
          if (temperature[i] < T_solidus)
            jacobian_vector[i] = 0.0;
          else if (temperature[i] > T_liquidus)
            jacobian_vector[i] = surface_tension_gradient;
          else
            {
              const double l_frac = calculate_liquid_fraction(temperature[i],
                                                              T_solidus,
                                                              T_liquidus);
              jacobian_vector[i]  = surface_tension_gradient * l_frac;
            }
      }
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double T_solidus;
  const double T_liquidus;
};

#endif
