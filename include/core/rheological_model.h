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

template <int dim>
class RheologicalModel
{
public: 

  RheologicalModel()
    {}

  virtual double get_viscosity(double shear_rate_magnitude) = 0;
  
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
  Carreau(Parameters::NonNewtonian non_newtonian_parameters)
    : viscosity_0(non_newtonian_parameters.viscosity_0)
    , viscosity_inf(non_newtonian_parameters.viscosity_inf)
    , lambda(non_newtonian_parameters.lambda)
    , a(non_newtonian_parameters.a)
    , n(non_newtonian_parameters.n)
    {}

  double get_viscosity(double shear_rate_magnitude) override
  {
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
