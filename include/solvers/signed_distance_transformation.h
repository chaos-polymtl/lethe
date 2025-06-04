// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_signed_distance_transformation_h
#define lethe_signed_distance_transformation_h


#include <core/parameters_multiphysics.h>

/**
 * @brief Abstract class for filtering phase fraction in volume of fluid (VOF)
 * simulations.
 */
class SignedDistanceTransformationBase
{
public:
  /**
   * @brief Instantiate and return a pointer to a SignedDistanceTransformationBase
   * object by casting it to the proper child class.
   *
   * @param[in] geometric_redistanciation_parameters VOF geometric redistanciation parameters.
   *
   * @return Casted SignedDistanceTransformationBase child class object.
   */
  static std::shared_ptr<SignedDistanceTransformationBase>
  model_cast(
    const Parameters::VOF_GeometricInterfaceReinitialization &geometric_redistanciation_parameters);

  /**
   * @brief Compute the value of the phase fraction from the signed distance.
   *
   * @param[in] signed_distance Value of the signed distance
   *
   * @return Value of the computed phase fraction after applying the transtion.
   */
  virtual double
  transfom_signed_distance(const double signed_distance) = 0;
};

/**
 * @brief Transforms the signed distance to the phase fraction with a hyperbolic tangent function.
 *
 * The  phase fraction is defined as:
 *
 * \f$\phi = 0.5 - 0.5 \tanh(d/\varepsilon)\f$
 *
 * where \f$d\f$ is the signed distance and \f$\varepsilon\f$
 * is a model parameter giving a measure of the interface thickness.
 *
 */
class SignedDistanceTransformationTanh : public SignedDistanceTransformationBase
{
public:
  /**
   * @brief Constructor of the hyperbolic tangent signed distance transformation.
   *
   * @param[in] tanh_thickness Value of the \f$\varepsilon\f$ parameter modulating the
   * interface thickness.
   */
  SignedDistanceTransformationTanh(const double tanh_thickness)
    : tanh_thickness(tanh_thickness)
  {}

  /**
   * @brief Computes the value of the phase fraction from the signed distance.
   *
   * @param[in] signed_distance Value of the signed distance
   *
   * @return Value of the computed phase fraction after applying the transtion.
   */
  double
  transfom_signed_distance(const double signed_distance) override
  {
    return 0.5 - 0.5 * tanh(signed_distance/tanh_thickness);
  }

private:
  /**
   * User-defined parameter that defines the thickness of the interface. It influences its sharpness: a smaller value leads to a sharper interface.
   */
  const double tanh_thickness;
};

/**
 * @brief Transforms the signed distance to the phase fraction with a piecewise polynomial function (degree 4).
 *
 * The  phase fraction is defined as:
 *
 * if \f$d < 0 \f$
 * \f$\phi = 0.5 - 0.5(4d' + 6d'^2 + 4*d'^3 + d'^4)\f$
 *
 * else
 * \f$\phi = 0.5 - 0.5(4d' - 6d'^2 + 4*d'^3 - d'^4)\f$
 * 
 * with \f$d' = d/d_\mathrm{max}\f$, where \f$d\f$ the signed distance and \f$d_\mathrm{max}\f$ the maximum redistanciation distance. This transformation clamps the phase fraction to 0 or 1 when \f$d = \pm d_\mathrm{max}\f$.
 *
 */
class SignedDistanceTransformationPiecewisePolynomial : public SignedDistanceTransformationBase
{
public:
  /**
   * @brief Constructor of the piecewise polynomial signed distance transformation.
   *
   * @param[in] max_reinitialization_distance Value of the maximum distance
   */
  SignedDistanceTransformationPiecewisePolynomial(const double max_reinitialization_distance)
    : max_reinitialization_distance(max_reinitialization_distance)
  {}

  /**
   * @brief Computes the value of the phase fraction from the signed distance.
   *
   * @param[in] signed_distance Value of the signed distance
   *
   * @return Value of the computed phase fraction after applying the transtion.
   */
  double
  transfom_signed_distance(const double signed_distance) override
  {
    const double dimentionless_d = signed_distance/max_reinitialization_distance;
    
    double piecewise_polynomial_value;
    
    if (dimentionless_d < 0.0)
    {
      piecewise_polynomial_value = 4.0*dimentionless_d + 6.0*Utilities::fixed_power<2>(dimentionless_d) + 4.0*Utilities::fixed_power<3>(dimentionless_d) + Utilities::fixed_power<4>(dimentionless_d);
    }
    else
    {
      piecewise_polynomial_value = 4.0*dimentionless_d - 6.0*Utilities::fixed_power<2>(dimentionless_d) + 4.0*Utilities::fixed_power<3>(dimentionless_d) - Utilities::fixed_power<4>(dimentionless_d);
    }
    return 0.5 - 0.5 * piecewise_polynomial_value;
  }

private:
  /**
   * Maximum redistanciation distance.
   */
  const double max_reinitialization_distance;
};

#endif
