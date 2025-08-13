// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/output_struct.h>

#include <deal.II/base/vector_tools>

template struct OutputStructPostprocessor<2, GlobalVectorType>;
template struct OutputStructPostprocessor<3, GlobalVectorType>;
template struct OutputStructPostprocessor<2, GlobalBlockVectorType>;
template struct OutputStructPostprocessor<3, GlobalBlockVectorType>;
template struct OutputStructPostprocessor<
  2,
  LinearAlgebra::distributed::Vector<double>>;
template struct OutputStructPostprocessor<
  3,
  LinearAlgebra::distributed::Vector<double>>;

template struct OutputStructSolution<2, GlobalVectorType>;
template struct OutputStructSolution<3, GlobalVectorType>;
template struct OutputStructSolution<2, GlobalBlockVectorType>;
template struct OutputStructSolution<3, GlobalBlockVectorType>;
template struct OutputStructSolution<
  2,
  LinearAlgebra::distributed::Vector<double>>;
template struct OutputStructDoFHandler<
  3,
  LinearAlgebra::distributed::Vector<double>>;
