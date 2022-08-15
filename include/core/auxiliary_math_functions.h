/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
#include "core/parameters.h"

#include <cfloat>

#ifndef LETHE_AUXILIARY_MATH_FUNCTIONS_H
#  define LETHE_AUXILIARY_MATH_FUNCTIONS_H

/**
 * Carries out the calculation of the harmonic mean of two values.
 * The harmonic mean of x and y is defined by 2*x*y/(x+y)
 *
 * @param value_one
 * @param value_two
 * @return harmonic mean of value_one and value_two
 */
inline double
harmonic_mean(const double &value_one, const double &value_two)
{
  return (2 * value_one * value_two / (value_one + value_two + DBL_MIN));
}

#endif // LETHE_AUXILIARY_MATH_FUNCTIONS_H
