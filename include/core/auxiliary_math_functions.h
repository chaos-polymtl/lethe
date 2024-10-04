// SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_auxiliary_math_functions_h
#define lethe_auxiliary_math_functions_h

#include <limits>

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
  return (2 * value_one * value_two /
          (value_one + value_two + std::numeric_limits<double>::min()));
}

#endif
