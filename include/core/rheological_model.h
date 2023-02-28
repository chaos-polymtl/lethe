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
   * @brief Instanciates and returns a pointer to a RheologicalModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parsed physical properties that will provide
   * either the model rheological model being used or say it is a
   * Newtonian flow
   */
  static std::shared_ptr<RheologicalModel>
  model_cast(const Parameters::Material &material_properties);

  /**
   * @brief Returns the value of the n parameters of the model, is the model has one.
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
   * @brief Returns the value of the viscosity, is the model has one.
   */
  virtual double
  get_viscosity() const
  {
    return 1.0;
  }

  /**
   * @brief Sets a new viscosity value, if the models has one.
   */
  virtual void
  set_viscosity(const double &)
  {}
};

class Newtonian : public RheologicalModel
{
public:
  /**
   * @brief Parameter constructor
   *
   * @param viscosity The constant newtonian viscosity
   */
  Newtonian(const double &p_viscosity)
    : viscosity(p_viscosity)
  {}

  /**
   * @brief Returns the viscosity.
   *
   * @param field_values Values of the field on which the viscosity may depend on.
   * For the constant viscosity, the viscosity does not depend on anything.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the viscosity.
   * @param field_vectors Values of the field on which the viscosity may depend on. These are not used for the constant viscosity
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(), property_vector.end(), viscosity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the viscosity
   * with respect to a field
   * @param field_values Value of the various fields on which the property may depend. The constant viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the viscosity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the viscosity with respect to the field id
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
  get_viscosity() const override
  {
    return viscosity;
  }

  void
  set_viscosity(const double &p_viscosity) override
  {
    viscosity = p_viscosity;
  }

private:
  double viscosity;
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
   * @brief Returns the non-newtonian viscosity.
   *
   * @param field_values The values Values of the field on which the viscosity may depend on. For this model, it only depends on the magnitude of the shear rate tensor
   *
   * Source : Morrison, F. A. (2001). No Memory: Generalized Newtonian Fluids.
   * Understanding Rheology. Raymond F. Boyer Librabry Collection, Oxford
   * University Press.
   */
  double
  value(const std::map<field, double> &field_values) override;

  /**
   * @brief vector_value Calculates the vector values of the viscosity.
   * @param field_vectors Values of the field on which the viscosity may depend on. The power-law viscosity depends on the shear rate.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the viscosity
   * with respect to a field
   * @param field_values Value of the various fields on which the property may depend. The constant viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the viscosity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the viscosity with respect to the field id
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
  get_viscosity() const override
  {
    return K;
  }

  void
  set_viscosity(const double &p_viscosity) override
  {
    K = p_viscosity;
  }

private:
  inline double
  calculate_viscosity(const double shear_rate_magnitude)
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
  Carreau(const double viscosity_0,
          const double viscosity_inf,
          const double lambda,
          const double a,
          const double n)
    : viscosity_0(viscosity_0)
    , viscosity_inf(viscosity_inf)
    , lambda(lambda)
    , a(a)
    , n(n)
  {
    this->model_depends_on[shear_rate] = true;
  }

  /**
   * @brief Returns the non-newtonian viscosity.
   *
   * @param field_values The values of the field on which the viscosity may depend on. For this model, it only depends on the magnitude of the shear rate tensor
   *
   * Source : Morrison, F. A. (2001). No Memory: Generalized Newtonian Fluids.
   * Understanding Rheology. Raymond F. Boyer Librabry Collection, Oxford
   * University Press.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the viscosity.
   * @param field_vectors Values of the field on which the viscosity may depend on. The power-law viscosity depends on the shear rate.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the viscosity
   * with respect to a field
   * @param field_values Value of the various fields on which the property may depend. The constant viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the viscosity with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the viscosity with respect to the field id
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
  get_viscosity() const override
  {
    return viscosity_0;
  }

  void
  set_viscosity(const double &p_viscosity) override
  {
    viscosity_0 = p_viscosity;
  }


private:
  inline double
  calculate_viscosity(const double shear_rate_magnitude)
  {
    return viscosity_inf +
           (viscosity_0 - viscosity_inf) *
             std::pow(1.0 + std::pow(shear_rate_magnitude * lambda, a),
                      (n - 1) / a);
  }

  double viscosity_0;
  double viscosity_inf;
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
   * @param p_phase_change parameters The parameters needed by the phase change rheology
   */
  PhaseChangeRheology(const Parameters::PhaseChange p_phase_change_parameters)
    : param(p_phase_change_parameters)
  {
    this->model_depends_on[temperature] = true;
  }

  /**
   * @brief Returns the viscosity.
   *
   * @param field_values Values of the field on which the viscosity may depend on.
   * For the constant viscosity, the viscosity does not depend on anything.
   */
  double
  value(const std::map<field, double> & /*field_values*/) override;

  /**
   * @brief vector_value Calculates the vector values of the viscosity.
   * @param field_vectors Values of the field on which the viscosity may depend on. These are not used for the constant viscosity
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override;

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the viscosity
   * with respect to a field
   * @param field_values Value of the various fields on which the property may depend. The constant viscosity does not depend on anything
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the viscosity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override;

  /**
   * @brief vector_jacobian Calculate the derivative of the with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the viscosity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override;

private:
  /**
   * @brief viscosity Calculates the viscosity from the temperature using the liquid fraction
   *
   * @param T value of the temperature
   * @return value of the viscosity
   *
   */
  inline double
  viscosity(const double T)
  {
    const double l_frac =
      std::min(std::max((T - param.T_solidus) /
                          (param.T_liquidus - param.T_solidus),
                        0.),
               1.);

    return param.viscosity_l * l_frac + param.viscosity_s * (1. - l_frac);
  }

  Parameters::PhaseChange param;
};

#endif
