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

using namespace dealii;

/**
 * @brief RheologicalModel. Abstract class that allows to calculate the 
 * non-newtonian viscosity on each quadrature point and the shear rate
 * magnitude. RheologicalModel::get_viscosity() is a pure virtual method,
 * since it can only be calculated knowing the rheological model that's
 * being used.
 */
template <int dim>
class RheologicalModel
{
public: 
  /**
   * @brief Default constructor
   */
  RheologicalModel()
    {}

  /**
   * @brief Returns the non-newtonian viscosity.
   *
   * @param shear_rate The shear rate tensor at the position of the 
   * considered quadratured point
   */
  virtual double get_viscosity(const Tensor<2, dim> shear_rate) = 0;
  
  /**
   * @brief Returns the magnitude of the shear rate tensor given in parameter.
   *
   * @param shear_rate The shear rate tensor at the position of the 
   * considered quadratured point
   */
  double 
  get_shear_rate_magnitude(const Tensor<2, dim> shear_rate)
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

  
};

template <int dim>
class Carreau : public RheologicalModel<dim>
{
public:
  /**
   * @brief Parameter constructor
   * 
   * @param non_newtonian_parameters The non newtonian parameters
   */
  Carreau(Parameters::NonNewtonian non_newtonian_parameters)
    : viscosity_0(non_newtonian_parameters.viscosity_0)
    , viscosity_inf(non_newtonian_parameters.viscosity_inf)
    , lambda(non_newtonian_parameters.lambda)
    , a(non_newtonian_parameters.a)
    , n(non_newtonian_parameters.n)
    {}

  /**
   * @brief Returns the non-newtonian viscosity.
   *
   * @param shear_rate The shear rate tensor at the position of the 
   * considered quadratured point
   */
  double 
  get_viscosity(const Tensor<2, dim> shear_rate) override
  {
    double shear_rate_magnitude = this->get_shear_rate_magnitude(shear_rate);
    return viscosity_inf + (viscosity_0 - viscosity_inf) * std::pow(1.0 + std::pow(shear_rate_magnitude * lambda, a), (n - 1)/a);
  }
  
private:
  double viscosity_0;
  double viscosity_inf;
  double lambda;
  double a;
  double n;
};


#endif
