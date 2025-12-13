// SPDX-FileCopyrightText: Copyright (c) 2021-2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
};


#endif // LETHE_RPT_CALCULATING_PARAMETERS_H
