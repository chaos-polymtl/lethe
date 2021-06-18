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
 *  This file defines parameters needed for the RPT simulation in the Parameters
 * namespace.
 */

#ifndef lethe_rpt_parameters_h
#define lethe_rpt_parameters_h

#include <deal.II/base/config.h>

#include <core/parameters.h>

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief RPTParameters - Defines the common parameters for the RPT
   * simulation as fixed parameters, number of Monte Carlo iterations and
   * filename of the particle positions.
   */

  struct RPTParameters
  {
    // Particle positions filename
    std::string particle_positions_file;

    // Enable to export counts result in a .csv file
    bool export_counts;

    // Exporting counts csv filename
    std::string export_counts_file;

    // Number of Monte Carlo iterations for alpha and theta
    unsigned int n_monte_carlo_iteration;

    // Seed of the random number generator
    int seed;

    // All parameters that are fixed by the user
    double reactor_radius; // [m]
    double peak_to_total_ratio;
    double sampling_time; // [s]

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief InitialRPTParameters - Allows parameters tuning. If turned on,
   * other parameters values are the initial guess for the optimization,
   * otherwise those act as fixed parameters.
   */

  struct InitialRPTParameters
  {
    // Enable tuning parameters
    bool tuning;

    // Filename of experimental data
    std::string experimental_file;

    // Parameters to tune or fixed parameters if tuning is disable
    double dead_time;          // Dead time of the detector per accepted pulse
    double activity;           // Activity of the tracer
    double gamma_rays_emitted; // Number of gamma-rays emitted by each
                               // disintegration
    double attenuation_coefficient_reactor;  // Total linear attenuation
                                             // coefficient of the medium
    double attenuation_coefficient_detector; // Total linear attenuation
                                             // coefficient of the detector

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief DetectorParameters - Defines information related to detectors. All
   * detectors must have the same dimensions (r and l).
   */

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
