// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_utilities_h
#define lethe_dem_utilities_h

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria_base.h>

using namespace dealii;


/**
 * @brief Reports the ratio between the minimum cell size and the maximum
 * particle diameter.
 *
 * This function loops over all locally owned cells in a distributed mesh
 * to find the minimum vertex-to-vertex distance (cell size), then compares
 * it to the maximum particle diameter. It prints out the results using
 * a ConditionalOStream and issues a warning if the ratio is too low.
 *
 * @tparam dim The spatial dimension of the triangulation.
 * @param triangulation A reference to the distributed triangulation (mesh) object.
 * @param maximum_particle_diameter The largest particle diameter in the simulation.
 * @param pcout A ConditionalOStream used for printing output. In parallel simulations,
 * this stream ensures that only one processor (usually the master) prints
 * messages to avoid cluttered output from all MPI ranks.
 * @param mpi_communicator The MPI communicator used for the parallel simulation.
 * This is required to compute the **global minimum cell size** across all
 * MPI processes.
 */
template <int dim>
void
report_cell_size_to_particle_diameter_ratio(
  const parallel::DistributedTriangulationBase<dim> &triangulation,
  const double                                       maximum_particle_diameter,
  const ConditionalOStream                          &pcout,
  const MPI_Comm                                    &mpi_communicator)
{
  double min_vertex_distance = std::numeric_limits<double>::max();
  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      min_vertex_distance =
        std::min(min_vertex_distance, cell->minimum_vertex_distance());
  min_vertex_distance =
    Utilities::MPI::min(min_vertex_distance, mpi_communicator);
  pcout << "Minimum vertex distance between cell vertices: "
        << min_vertex_distance << std::endl;
  double ratio = min_vertex_distance / maximum_particle_diameter;
  pcout << "Minimum cell size to maximum particle diameter ratio: " << ratio
        << std::endl;

  if (ratio < 1.0)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Minimum cell size is smaller than the maximum particle diameter. "
          "Consider coarsening the mesh to achieve a ratio larger than 1"));
    }
}


#endif
