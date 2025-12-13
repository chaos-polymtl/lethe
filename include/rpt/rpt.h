// SPDX-FileCopyrightText: Copyright (c) 2021-2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_rpt_h
#define lethe_rpt_h

/**
 * This class handles positions of particles and detector and allows to
 * calculate count for every pair. It also shows results in terminal or exports
 * them in .csv file if enable.
 */

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
   * @brief Set up data assignation & call photon counts calculation
   * for all particle positions and detector
   */
  void
  setup_and_calculate();

  /**
   * @brief Export position, detector id and count results in .csv or .dat
   */
  void
  export_data();

  /**
   * @brief Calculate photon counts for all particle positions and detector
   */
  void
  calculate_counts();

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
  std::vector<double>             calculated_counts;
};

#endif
