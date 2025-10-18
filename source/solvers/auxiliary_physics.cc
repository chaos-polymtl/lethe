// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/auxiliary_physics.h>

template <int dim, typename VectorType>
void
AuxiliaryPhysics<dim, VectorType>::serialize_tables_vector(
  const std::vector<OutputStructTableHandler> &table_output_structs,
  MPI_Comm                                     mpi_communicator)
{
  for (const auto &output_table : table_output_structs)
    {
      serialize_table(output_table.table,
                      output_table.table_filename,
                      mpi_communicator);
    }
}

template <int dim, typename VectorType>
void
AuxiliaryPhysics<dim, VectorType>::deserialize_tables_vector(
  std::vector<OutputStructTableHandler> &table_output_structs,
  MPI_Comm                               mpi_communicator)
{
  for (auto &output_table : table_output_structs)
    {
      deserialize_table(output_table.table,
                        output_table.table_filename,
                        mpi_communicator);
    }
}
// Pre-compile the 2D and 3D version with the types that can occur
template class AuxiliaryPhysics<2, GlobalVectorType>;
template class AuxiliaryPhysics<3, GlobalVectorType>;
template class AuxiliaryPhysics<2, GlobalBlockVectorType>;
template class AuxiliaryPhysics<3, GlobalBlockVectorType>;

#ifndef LETHE_USE_LDV
template class AuxiliaryPhysics<2, LinearAlgebra::distributed::Vector<double>>;
template class AuxiliaryPhysics<3, LinearAlgebra::distributed::Vector<double>>;
#endif
