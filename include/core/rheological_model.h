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

#ifndef lethe_rheological_model_h
#define lethe_rheological_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief RheologicalModel. Abstract class that allows to calculate the
 * non-newtonian viscosity on each quadrature point and the shear rate
 * magnitude.
 */
class RheologicalModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Default constructor
   */
  RheologicalModel()
  {}

  /**
   * @brief Instanciates and returns a pointer to a RheologicalModel object by
   * casting it to
   * the proper child class
   *
   * @param material_properties Parsed physical properties that will provide
   * either the model rheological model being used or say it is a
   * Newtonian flow
   */
  static std::shared_ptr<RheologicalModel>
  model_cast(const Parameters::Material &material_properties);

  /**
   * @brief Returns the value of the n parameters of the model, if the model
   * has one.
   */
  virtual double
  get_n() const
  {
    return 1.0;
  }

  /**
   * @brief Sets a new value to the n parameter, if the models has one.
   */
  virtual void
  set_n(const double &)
  {}

  /**
   * @brief Returns the value of the kinematic viscosity, if the model has one.
   */
  virtual double
  get_kinematic_viscosity() const
  {
    return 1.0;
  }

  /**
   * @brief Returns the kinematic viscosity scale. This is used to get the
   * kinematic viscosity scale in the problem when it's not constant. This is
   * useful for weak boundary conditions.
   */
  virtual double
  get_kinematic_viscosity_scale() const
  {
    return 1.0;
  }

  /**
   * @brief Sets a new kinematic viscosity value, if the models has one.
   * @param p_kinematic_viscosity The kinematic viscosity value
   */
  virtual void
  set_kinematic_viscosity(const double & /*p_kinematic_viscosity*/)
  {}

  /**
   * @brief Returns the value of the dynamic viscosity, if a kinematic viscosity
   * and a reference density value is specified.
   * @param p_density_ref The density of the fluid at the reference state
   * @param p_shear_rate_magnitude The shear rate magnitude
   * @param p_temperature The temperature of the fluid
   */
  virtual double
  get_dynamic_viscosity(const double & /*p_density_ref*/,
                        const double & /*p_shear_rate_magnitude*/,
                        const double & /*p_temperature*/) const = 0;
};

class Newtonian : public RheologicalModel
{
public:
  /**
   * @brief Parameter constructor
   *
   * @param p_kinematic_viscosity The constant newtonian viscosity
   */
  Newtonian(const double &p_kinematic_viscosity)
    : kinematic_viscosity(p_kinematic_viscosity)
  {}

  /**
   * @brief Returns the kinematic viscosity.
   *
   * @param field_values Values of the field on which the viscosity may depend
   * on.
   * For the constant kinematic viscosity, the kinematic viscosity does not
   * depend on anything.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the kinematic
   * viscosity.
   * @param field_vectors Values of the field on which the kinematic viscosity
   * may depend on. These are not used for the constant kinematic viscosity
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              kinematic_viscosity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * kinematic viscosity
   * with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend. The constant kinematic viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the kinematic viscosity with
   * respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a
   * field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the kinematic
   * viscosity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

  double
  get_kinematic_viscosity() const override
  {
    return kinematic_viscosity;
  }

  double
  get_kinematic_viscosity_scale() const override
  {
    return kinematic_viscosity;
  }

  void
  set_kinematic_viscosity(const double &p_kinematic_viscosity) override
  {
    kinematic_viscosity = p_kinematic_viscosity;
  }

  double
  get_dynamic_viscosity(const double &p_density_ref,
                        const double & /*p_shear_rate_magnitude*/,
                        const double & /*p_temperature*/) const override
  {
    return kinematic_viscosity / p_density_ref;
  }

private:
  double kinematic_viscosity;
};

class PowerLaw : public RheologicalModel
{
public:
  /**
   * @brief Parameter constructor
   *
   * @param non_newtonian_parameters The non newtonian parameters
   */
  PowerLaw(const double K, const double n, const double shear_rate_min)
    : K(K)
    , n(n)
    , shear_rate_min(shear_rate_min)
  {
    this->model_depends_on[shear_rate] = true;
  }

  /**
   * @brief Returns the non-newtonian kinematic viscosity.
   *
   * @param field_values The values Values of the field on which the kinematic
   * viscosity may depend on. For this model, it only depends on the magnitude
   * of the shear rate tensor
   *
   * Source : Morrison, F. A. (2001). No Memory: Generalized Newtonian Fluids.
   * Understanding Rheology. Raymond F. Boyer Librabry Collection, Oxford
   * University Press.
   */
  double
  value(const std::map<field, double> &field_values) override;

  /**
   * @brief vector_value Calculates the vector values of the kinematic
   * viscosity.
   * @param field_vectors Values of the field on which the kinematic viscosity
   * may depend on. The power-law kinematic viscosity depends on the shear rate.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * kinematic viscosity with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend. The constant kinematic viscosity does not depend on anything.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the kinematic viscosity with
   * respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a
   * field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property
   * @param id Identifier of the field with respect to which a derivative
   * should be calculated
   * @param jacobian Vector of the value of the derivative of the kinematic
   * viscosity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override;

  double
  get_n() const override
  {
    return n;
  }

  void
  set_n(const double &p_n) override
  {
    n = p_n;
  }

  double
  get_kinematic_viscosity() const override
  {
    return K;
  }

  double
  get_kinematic_viscosity_scale() const override
  {
    return 1.0 > shear_rate_min ? K * std::pow(1.0, n - 1) :
                                  K * std::pow(shear_rate_min, n - 1);
  }

  void
  set_kinematic_viscosity(const double &p_kinematic_viscosity) override
  {
    K = p_kinematic_viscosity;
  }

  double
  get_dynamic_viscosity(const double &p_density_ref,
                        const double &p_shear_rate_magnitude,
                        const double & /*p_temperature*/) const override
  {
    return calculate_kinematic_viscosity(p_shear_rate_magnitude) /
           p_density_ref;
  }

private:
  inline double
  calculate_kinematic_viscosity(const double shear_rate_magnitude) const
  {
    return shear_rate_magnitude > shear_rate_min ?
             K * std::pow(shear_rate_magnitude, n - 1) :
             K * std::pow(shear_rate_min, n - 1);
  }

  inline double
  calculate_derivative(const double shear_rate_magnitude)
  {
    return shear_rate_magnitude > shear_rate_min ?
             (n - 1) * K * std::pow(shear_rate_magnitude, n - 2) :
             0;
  }

  double       K;
  double       n;
  const double shear_rate_min;
};

class Carreau : public RheologicalModel
{
public:
  /**
   * @brief Parameter constructor
   *
   * @param non_newtonian_parameters The non newtonian parameters
   */
  Carreau(const double p_kinematic_viscosity_0,
          const double p_kinematicviscosity_inf,
          const double p_lambda,
          const double p_a,
          const double p_n)
    : kinematic_viscosity_0(p_kinematic_viscosity_0)
    , kinematic_viscosity_inf(p_kinematicviscosity_inf)
    , lambda(p_lambda)
    , a(p_a)
    , n(p_n)
  {
    this->model_depends_on[shear_rate] = true;
  }

  /**
   * @brief Returns the non-newtonian viscosity.
   *
   * @param field_values The values of the field on which the viscosity may
   * depend on. For this model, it only depends on the magnitude of the shear
   * rate tensor.
   *
   * Source : Morrison, F. A. (2001). No Memory: Generalized Newtonian Fluids.
   * Understanding Rheology. Raymond F. Boyer Librabry Collection, Oxford
   * University Press.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the viscosity.
   * @param field_vectors Values of the field on which the viscosity may depend
   * on. The power-law viscosity depends on the shear rate.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * viscosity with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend. The constant viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the viscosity with respect to
   * the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a
   * field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the viscosity with
   * respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override;

  double
  get_n() const override
  {
    return n;
  }

  void
  set_n(const double &p_n) override
  {
    n = p_n;
  }

  double
  get_kinematic_viscosity() const override
  {
    return kinematic_viscosity_0;
  }

  double
  get_kinematic_viscosity_scale() const override
  {
    return kinematic_viscosity_inf +
           (kinematic_viscosity_0 - kinematic_viscosity_inf) *
             std::pow(1.0 + std::pow(1.0 * lambda, a), (n - 1) / a);
  }

  void
  set_kinematic_viscosity(const double &p_kinematic_viscosity) override
  {
    kinematic_viscosity_0 = p_kinematic_viscosity;
  }

  double
  get_dynamic_viscosity(const double &p_density_ref,
                        const double &p_shear_rate_magnitude,
                        const double & /*p_temperature*/) const override
  {
    return calculate_kinematic_viscosity(p_shear_rate_magnitude) /
           p_density_ref;
  }

private:
  inline double
  calculate_kinematic_viscosity(const double shear_rate_magnitude) const
  {
    return kinematic_viscosity_inf +
           (kinematic_viscosity_0 - kinematic_viscosity_inf) *
             std::pow(1.0 + std::pow(shear_rate_magnitude * lambda, a),
                      (n - 1) / a);
  }

  double kinematic_viscosity_0;
  double kinematic_viscosity_inf;
  double lambda;
  double a;
  double n;
};

/**
 * @brief Calculates the magnitude of the shear rate tensor given as input
 *
 * @param shear_rate A shear rate tensor
 */
template <int dim>
inline double
calculate_shear_rate_magnitude(const Tensor<2, dim> shear_rate)
{
  double shear_rate_magnitude = 0;
  for (unsigned int i = 0; i < dim; ++i)
    {
      for (unsigned int j = 0; j < dim; ++j)
        {
          shear_rate_magnitude += (shear_rate[i][j] * shear_rate[j][i]);
        }
    }
  shear_rate_magnitude = sqrt(0.5 * shear_rate_magnitude);
  return shear_rate_magnitude;
}


class PhaseChangeRheology : public RheologicalModel
{
public:
  /**
   * @brief Parameter constructor
   *
   * @param p_phase_change parameters The parameters needed by the phase change
   * rheology
   */
  PhaseChangeRheology(const Parameters::PhaseChange p_phase_change_parameters)
    : param(p_phase_change_parameters)
  {
    this->model_depends_on[temperature] = true;
  }

  /**
   * @brief Returns the kinematic viscosity.
   *
   * @param field_values Values of the field on which the kinematic viscosity
   * may depend on.
   * For the constant kinematic viscosity, the kinematic viscosity does not
   * depend on anything.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the kinematic viscosity.
   * @param field_vectors Values of the field on which the kinematic viscosity
   * may depend on. These are not used for the constant kinematic viscosity.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the
   * kinematic viscosity with respect to a field
   * @param field_values Value of the various fields on which the property may
   * depend. The constant kinematic viscosity does not depend on anything.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the kinematic viscosity with
   * respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a
   * field
   * @param field_vectors Vector for the values of the fields used to evaluate
   * the property
   * @param id Identifier of the field with respect to which a derivative should
   * be calculated
   * @param jacobian Vector of the value of the derivative of the kinematic
   * viscosity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override;

  double
  get_dynamic_viscosity(const double &p_density_ref,
                        const double & /*p_shear_rate_magnitude*/,
                        const double &p_temperature) const override
  {
    return kinematic_viscosity(p_temperature) / p_density_ref;
  }

private:
  /**
   * @brief kinematic_viscosity Calculates the kinematic viscosity from the
   * temperature using the liquid fraction
   *
   * @param T value of the temperature
   * @return value of the kinematic viscosity
   *
   */
  inline double
  kinematic_viscosity(const double T) const
  {
    const double l_frac =
      std::min(std::max((T - param.T_solidus) /
                          (param.T_liquidus - param.T_solidus),
                        0.),
               1.);

    return param.kinematic_viscosity_l * l_frac + param.kinematic_viscosity_s * (1. - l_frac);
  }

  Parameters::PhaseChange param;
};

#endif
