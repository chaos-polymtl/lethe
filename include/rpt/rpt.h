/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Authors: Bruno Blais, Ghazaleh Mirakhori, Audrey Collard-Daigneault,
* Polytechnique Montreal, 2020-
*/

/**
 * This class allows to calculate the photon count from a particle received by
 * a detector with the Monte Carlo method.
 */

#ifndef lethe_rpt_h
#define lethe_rpt_h

/**
 * This class handles positions of particles and detector and allows to
 * calculate count for every pair. It also shows results in terminal or exports
 * them in .csv file if enable.
 */

// deal.II includes
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

using namespace dealii;

template <int dim>
class RPT
{
public:
  /**
   * @brief Constructor for the RPT
   *
   * @param RPTparameters All parameters and positions needed for the count
   * calculation
   *
   */
  RPT(RPTCalculatingParameters &RPTparameters);

  /**
   * @brief Calculate photon count for all particle positions and detector
   */
  void
  calculate();

  /**
   * @brief Read text file for particle positions and assign them to particle
   * objects
   */
  void
  assign_particle_positions();

  /**
   * @brief Read text file for detector positions and assign them to detector
   * objects
   */
  void
  assign_detector_positions();

  /**
   * @brief Read text file for experimental count data and store it in a vector
   */
  std::vector<double>
  extract_experimental_counts();

  /**
   * @brief Calculate cost function for calculated and measured counts
   */
  double
  calculate_cost_function(std::vector<double> &measured_counts,
                          std::vector<double> &calculated_counts);

  RPTCalculatingParameters rpt_parameters;

private:
  std::vector<RadioParticle<dim>> particle_positions;
  std::vector<Detector<dim>>      detectors;
};

#endif
