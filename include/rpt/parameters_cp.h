/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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

/**
 *  This file defines parameters needed for the CP simulation in the Parameters
 * namespace.
 */

#ifndef lethe_cp_parameters_h
#define lethe_cp_parameters_h

#include <deal.II/base/config.h>

#include <core/parameters.h>

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief CPParameters - Defines the common parameters for the CP
   */

  struct CPParameters
  {
    int subdivisions;
    double cylinder_radius;
    double cylinder_half_length;
    unsigned int initial_refinement;
    double electric_field_tolerance;
    double velocity_tolerance;
    double CFD_input_velocity;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

} // namespace Parameters

class CPCalculatingParameters
{
public:
    CPCalculatingParameters()
    {}

    void
    declare(ParameterHandler &prm);

    void 
    parse(ParameterHandler &prm);

    Parameters::CPParameters cp_param;
};

#endif // lethe_cp_parameters_h