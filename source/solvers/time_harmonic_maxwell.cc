// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell.h>

using VectorType = GlobalVectorType;

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_dofs()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::percolate_time_vectors()
{
  // No time-dependent vectors to percolate in time-harmonic Maxwell
}


template class TimeHarmonicMaxwell<2>;
template class TimeHarmonicMaxwell<3>;
