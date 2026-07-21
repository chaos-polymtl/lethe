// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vector_h
#define lethe_vector_h

#include <deal.II/base/mpi.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector_operation.h>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>

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
 * @brief Traits describing a distributed vector type.
 *
 * These allow code that is generic with respect to the vector type to branch
 * on whether a vector is a block vector, and to name the type of an individual
 * block, without having to know whether Lethe was configured with Trilinos
 * vectors or with deal.II distributed vectors.
 *
 * @tparam VectorType Type of the distributed vector.
 */
template <typename VectorType>
struct LetheVectorTraits
{
  /**
   * @brief Whether the vector is a deal.II distributed vector.
   *
   * This must never be expressed as a comparison against GlobalVectorType:
   * when Lethe is configured with LETHE_USE_LDV, GlobalVectorType *is* a
   * deal.II distributed vector, and such a comparison silently selects the
   * Trilinos behavior.
   */
  static constexpr bool is_dealii_vector =
    std::is_same_v<VectorType,
                   dealii::LinearAlgebra::distributed::Vector<double>> ||
    std::is_same_v<VectorType,
                   dealii::LinearAlgebra::distributed::BlockVector<double>>;

  static constexpr bool is_block_vector =
    std::is_same_v<VectorType, dealii::TrilinosWrappers::MPI::BlockVector> ||
    std::is_same_v<VectorType,
                   dealii::LinearAlgebra::distributed::BlockVector<double>>;

  using UnderlyingVectorType = VectorType;
};

template <>
struct LetheVectorTraits<dealii::TrilinosWrappers::MPI::BlockVector>
{
  static constexpr bool is_dealii_vector = false;
  static constexpr bool is_block_vector  = true;
  using UnderlyingVectorType             = dealii::TrilinosWrappers::MPI::Vector;
};

template <>
struct LetheVectorTraits<
  dealii::LinearAlgebra::distributed::BlockVector<double>>
{
  static constexpr bool is_dealii_vector = true;
  static constexpr bool is_block_vector  = true;
  using UnderlyingVectorType =
    dealii::LinearAlgebra::distributed::Vector<double>;
};

/**
 * @brief Reinitialize a vector that is written into during assembly.
 *
 * AffineConstraints::distribute_local_to_global and VectorTools::interpolate
 * write to degrees of freedom that are not owned by the current process.
 * Trilinos vectors buffer such off-processor writes internally and must
 * therefore stay free of ghost entries, since a ghosted Trilinos vector is
 * read-only. deal.II distributed vectors have no such buffer and can only
 * accept those writes if the corresponding entries are allocated as ghost
 * entries; the subsequent compress(VectorOperation::add) sends them to their
 * owner and clears them.
 *
 * Note that the ghost set has to cover more than the degrees of freedom of the
 * locally owned cells: a contribution to a constrained degree of freedom is
 * redistributed onto the degrees of freedom that constrain it, which for a
 * periodic boundary lies on the opposite side of the domain. Passing the set
 * the constraints were built with is the safe choice.
 *
 * @tparam VectorType Type of the distributed vector.
 * @tparam DofsType Type of storage of the indices (IndexSet for regular
 * vectors, std::vector<IndexSet> for block vectors).
 * @param[out] vector Vector to reinitialize.
 * @param locally_owned_dofs Degrees of freedom owned by this process.
 * @param ghost_dofs Degrees of freedom that may be written to on this process.
 * @param mpi_communicator MPI communicator of the vector.
 */
template <typename VectorType, typename DofsType>
void
reinit_assembly_vector(VectorType     &vector,
                       const DofsType &locally_owned_dofs,
                       const DofsType &ghost_dofs,
                       const MPI_Comm  mpi_communicator)
{
  if constexpr (LetheVectorTraits<VectorType>::is_dealii_vector)
    vector.reinit(locally_owned_dofs, ghost_dofs, mpi_communicator);
  else
    {
      (void)ghost_dofs;
      vector.reinit(locally_owned_dofs, mpi_communicator);
    }
}

/**
 * @brief Reinitialize a ghosted vector that is read during assembly.
 *
 * A freshly reinitialized deal.II distributed vector has its ghost entries
 * allocated but not imported, and reading one of them is an error until
 * update_ghost_values() has been called at least once. Trilinos vectors have
 * no such state and their update_ghost_values() is an empty function, so this
 * is free for them.
 *
 * @tparam VectorType Type of the distributed vector.
 * @tparam DofsType Type of storage of the indices (IndexSet for regular
 * vectors, std::vector<IndexSet> for block vectors).
 * @param[out] vector Vector to reinitialize.
 * @param locally_owned_dofs Degrees of freedom owned by this process.
 * @param locally_relevant_dofs Degrees of freedom readable by this process.
 * @param mpi_communicator MPI communicator of the vector.
 */
template <typename VectorType, typename DofsType>
void
reinit_ghosted_vector(VectorType     &vector,
                      const DofsType &locally_owned_dofs,
                      const DofsType &locally_relevant_dofs,
                      const MPI_Comm  mpi_communicator)
{
  vector.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  vector.update_ghost_values();
}

/**
 * @brief Get the range of indices locally owned by the current MPI process.
 *
 * Trilinos vectors and deal.II distributed vectors expose this information
 * through different interfaces, hence this compatibility helper.
 *
 * @tparam VectorType Type of the distributed vector.
 * @param vector Vector whose locally owned range is queried.
 * @return Half-open range [first, second) of locally owned indices.
 */
template <typename VectorType>
std::pair<dealii::types::global_dof_index, dealii::types::global_dof_index>
get_local_range(const VectorType &vector)
{
  if constexpr (std::is_same_v<VectorType,
                               dealii::TrilinosWrappers::MPI::Vector>)
    return vector.local_range();
  else
    return vector.get_partitioner()->local_range();
}

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

/**
 * @brief Compute the global maximum entry of a distributed vector.
 *
 * Both Trilinos vectors and deal.II distributed vectors are supported. This
 * helper is required because dealii::LinearAlgebra::distributed::Vector does
 * not provide the max() and min() methods of TrilinosWrappers::MPI::Vector.
 *
 * @tparam VectorType Type of the distributed vector.
 * @param vector Vector over which the maximum is computed.
 * @param mpi_communicator MPI communicator over which the reduction is done.
 */
template <typename VectorType>
double
global_max(const VectorType &vector, const MPI_Comm mpi_communicator)
{
  double local_max = std::numeric_limits<double>::lowest();
  for (const auto i : vector.locally_owned_elements())
    local_max = std::max(local_max, vector[i]);

  return dealii::Utilities::MPI::max(local_max, mpi_communicator);
}

/**
 * @brief Compute the global minimum entry of a distributed vector.
 *
 * Both Trilinos vectors and deal.II distributed vectors are supported. This
 * helper is required because dealii::LinearAlgebra::distributed::Vector does
 * not provide the max() and min() methods of TrilinosWrappers::MPI::Vector.
 *
 * @tparam VectorType Type of the distributed vector.
 * @param vector Vector over which the minimum is computed.
 * @param mpi_communicator MPI communicator over which the reduction is done.
 */
template <typename VectorType>
double
global_min(const VectorType &vector, const MPI_Comm mpi_communicator)
{
  double local_min = std::numeric_limits<double>::max();
  for (const auto i : vector.locally_owned_elements())
    local_min = std::min(local_min, vector[i]);

  return dealii::Utilities::MPI::min(local_min, mpi_communicator);
}

#endif
