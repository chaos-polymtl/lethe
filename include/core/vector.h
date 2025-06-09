// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vector_h
#define lethe_vector_h

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operation.h>

// #define LETHE_USE_LDV is defined using the CMAKE flag -DLETHE_USE_LDV=TRUE

#ifndef LETHE_USE_LDV
using GlobalVectorType      = dealii::TrilinosWrappers::MPI::Vector;
using GlobalBlockVectorType = dealii::TrilinosWrappers::MPI::BlockVector;
#else
using GlobalVectorType = dealii::LinearAlgebra::distributed::Vector<double>;
using GlobalBlockVectorType =
  dealii::LinearAlgebra::distributed::BlockVector<double>;
#endif


/**
 * @brief Helper function that allows to convert deal.II vectors to Trilinos vectors.
 *
 * @tparam number Abstract type for number across the class (i.e., double).
 * @param out Destination TrilinosWrappers::MPI::Vector vector.
 * @param in Source LinearAlgebra::distributed::Vector<number> vector.
 */
template <typename number>
void
convert_vector_dealii_to_trilinos(
  dealii::TrilinosWrappers::MPI::Vector                    &out,
  const dealii::LinearAlgebra::distributed::Vector<number> &in)
{
  dealii::LinearAlgebra::ReadWriteVector<double> rwv(
    out.locally_owned_elements());
  rwv.import_elements(in, dealii::VectorOperation::insert);
  out.import_elements(rwv, dealii::VectorOperation::insert);
}

/**
 * @brief Helper function that allows to convert Trilinos vectors to deal.II vectors.
 *
 * @tparam number Abstract type for number across the class (i.e., double).
 * @param out Destination LinearAlgebra::distributed::Vector<number> vector.
 * @param in Source TrilinosWrappers::MPI::Vector vector.
 */
template <typename number>
void
convert_vector_trilinos_to_dealii(
  dealii::LinearAlgebra::distributed::Vector<number> &out,
  const dealii::TrilinosWrappers::MPI::Vector        &in)
{
  dealii::LinearAlgebra::ReadWriteVector<double> rwv(
    out.locally_owned_elements());
  rwv.reinit(in);
  out.import_elements(rwv, dealii::VectorOperation::insert);
}

#endif
