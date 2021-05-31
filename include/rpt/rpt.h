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
* Authors: Bruno Blais, Ghazaleh Mirakhori, Audrey Collard-Daigneault
Polytechnique Montreal, 2020-
*/


#ifndef lethe_rpt_h
#define lethe_rpt_h

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
  RPT(RPTCalculatingParameters &RPTparameters);

  void
  calculate();

  void
  assign_particle_positions();

  void
  assign_detector_positions();

private:
  RPTCalculatingParameters rpt_parameters;

  std::vector<RadioParticle<dim>> particle_positions;
  std::vector<Detector<dim>>      detectors;
};

#endif
