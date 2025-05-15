// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vector_h
#define lethe_vector_h

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// #define LETHE_USE_LDV is defined using the CMAKE flag -DLETHE_USE_LDV=TRUE

#ifndef LETHE_USE_LDV
using GlobalVectorType      = dealii::TrilinosWrappers::MPI::Vector;
using GlobalBlockVectorType = dealii::TrilinosWrappers::MPI::BlockVector;
#else
using GlobalVectorType = dealii::LinearAlgebra::distributed::Vector<double>;
using GlobalBlockVectorType =
  dealii::LinearAlgebra::distributed::BlockVector<double>;
#endif

#endif
