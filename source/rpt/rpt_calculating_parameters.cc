// SPDX-FileCopyrightText: Copyright (c) 2021-2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <rpt/rpt_calculating_parameters.h>

void
RPTCalculatingParameters::declare(ParameterHandler &prm)
{
  prm.declare_entry("print parameters",
                    "none",
                    Patterns::Selection("none|only changed|all"),
                    "Print all the parameters, or only"
                    "the changed parameters or none");

  Parameters::RPTParameters::declare_parameters(prm);
  Parameters::RPTTuningParameters::declare_parameters(prm);
  Parameters::DetectorParameters::declare_parameters(prm);
  Parameters::RPTReconstructionParameters::declare_parameters(prm);
  Parameters::RPTFEMReconstructionParameters::declare_parameters(prm);
}

void
RPTCalculatingParameters::parse(ParameterHandler &prm)
{
  rpt_param.parse_parameters(prm);
  tuning_param.parse_parameters(prm);
  detector_param.parse_parameters(prm);
  reconstruction_param.parse_parameters(prm);
  fem_reconstruction_param.parse_parameters(prm);
}
