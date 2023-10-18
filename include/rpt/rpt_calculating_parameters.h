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

#ifndef lethe_rpt_calculating_parameters_h
#define lethe_rpt_calculating_parameters_h

/**
 * Handles all the parameters declared in the parameter handler file for the
 * RPT.
 */

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>

#include <rpt/parameters_rpt.h>

class RPTCalculatingParameters
{
public:
  RPTCalculatingParameters()
  {}

  void
  declare(ParameterHandler &prm);

  void
  parse(ParameterHandler &prm);

  Parameters::RPTParameters                  rpt_param;
  Parameters::RPTTuningParameters            tuning_param;
  Parameters::DetectorParameters             detector_param;
  Parameters::RPTReconstructionParameters    reconstruction_param;
  Parameters::RPTFEMReconstructionParameters fem_reconstruction_param;
  Parameters::Mesh                           mesh;
};


#endif // LETHE_RPT_CALCULATING_PARAMETERS_H
