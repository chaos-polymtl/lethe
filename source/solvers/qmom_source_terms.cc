// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters.h>
#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

#include <span>

void
apply_qmom_sources(
  [[maybe_unused]] std::span<double>                      state,
  [[maybe_unused]] const double                           time_step,
  [[maybe_unused]] const Parameters::QMOM                &qmom_parameters,
  [[maybe_unused]] const Parameters::QMOMPhaseProperties &phase_properties)
{
  // The phase-space sources (nucleation, growth and coalescence) are not
  // implemented yet: the state is left unchanged.
}
