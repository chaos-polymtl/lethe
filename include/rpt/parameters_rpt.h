/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2021 by the Lethe authors
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

 *
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2021-
 */

/*
 * ---
 */

#ifndef lethe_rpt_parameters_h
#define lethe_rpt_parameters_h

#include <deal.II/base/config.h>

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <core/parameters.h>

using namespace dealii;

namespace Parameters
{
  struct RPTParameters
  {
    // File name with particle positions
    std::string particle_positions_file;

    // Number of Monte Carlo iteration for alpha and theta
    unsigned int iteration_number;

    // All parameters that are fixed by the user
    double reactor_radius;
    double peak_to_total_ratio;
    double sampling_time;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct InitialRPTParameters
  {
    bool   tuning;
    double dead_time;
    double activity;
    double gamma_rays_emitted;
    double attenuation_coefficient_reactor;
    double attenuation_coefficient_detector;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  // Detector parameters, all detector should have the same dimensions(r and l)
  struct DetectorParameters
  {
    double      radius;
    double      length;
    std::string detector_positions_file;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };



} // namespace Parameters



#endif // lethe_rpt_parameters_h
