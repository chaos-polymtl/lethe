// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/output_struct.h>
#include <core/vector.h>

#include <deal.II/lac/la_parallel_vector.h>

template struct OutputStructPostprocessor<2, GlobalVectorType>;
template struct OutputStructPostprocessor<3, GlobalVectorType>;
template struct OutputStructPostprocessor<2, GlobalBlockVectorType>;
template struct OutputStructPostprocessor<3, GlobalBlockVectorType>;
#ifndef LETHE_USE_LDV
template struct OutputStructPostprocessor<
  2,
  LinearAlgebra::distributed::Vector<double>>;
template struct OutputStructPostprocessor<
  3,
  LinearAlgebra::distributed::Vector<double>>;
#endif

template struct OutputStructSolution<2, GlobalVectorType>;
template struct OutputStructSolution<3, GlobalVectorType>;
template struct OutputStructSolution<2, GlobalBlockVectorType>;
template struct OutputStructSolution<3, GlobalBlockVectorType>;

#ifndef LETHE_USE_LDV
template struct OutputStructSolution<
  2,
  LinearAlgebra::distributed::Vector<double>>;
template struct OutputStructSolution<
  3,
  LinearAlgebra::distributed::Vector<double>>;
#endif
