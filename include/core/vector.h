/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 -  by the Lethe authors
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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_vector_h
#define lethe_vector_h

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#if false
using GlobalVectorType = dealii::TrilinosWrappers::MPI::Vector;
using GlobalBlockVectorType = dealii::TrilinosWrappers::MPI::BlockVector;
#else
using GlobalVectorType = dealii::LinearAlgebra::distributed::Vector<double>;
using GlobalBlockVectorType =
  dealii::LinearAlgebra::distributed::BlockVector<double>;
#endif

#endif
