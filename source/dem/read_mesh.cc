// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>

#include <dem/read_mesh.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <numbers>

template <int dim, int spacedim>
void
read_mesh(const Parameters::Mesh   &mesh_parameters,
          const bool                restart,
          const ConditionalOStream &pcout,
          parallel::DistributedTriangulationBase<dim, spacedim> &triangulation,
          const Parameters::Lagrangian::BCDEM                   &bc_params)
{
  pcout << "Reading triangulation" << std::endl;

  attach_grid_to_triangulation(triangulation, mesh_parameters);

  // Check if there are periodic boundaries
  for (unsigned int i_bc = 0; i_bc < bc_params.bc_types.size(); ++i_bc)
    {
      if (bc_params.bc_types[i_bc] ==
          Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
        match_periodic_boundaries(triangulation, bc_params);
    }

  if (!restart)
    {
      if (mesh_parameters.refine_until_target_size)
        {
          double minimal_cell_size =
            GridTools::minimal_cell_diameter(triangulation);
          double       target_size       = mesh_parameters.target_size;
          unsigned int number_refinement = static_cast<unsigned int>(std::floor(
            std::log(minimal_cell_size / target_size) / std::numbers::ln2));
          pcout << "Automatically refining grid until target size: "
                << target_size << std::endl;
          triangulation.refine_global(number_refinement);
          pcout << "Mesh was automatically refined: " << number_refinement
                << " times" << std::endl;
        }
      else
        {
          const int initial_refinement = mesh_parameters.initial_refinement;
          triangulation.refine_global(initial_refinement);
        }
    }

  pcout << std::endl << "Finished reading triangulation" << std::endl;
}

template <int dim, int spacedim>
void
match_periodic_boundaries(Triangulation<dim, spacedim>        &triangulation,
                          const Parameters::Lagrangian::BCDEM &bc_param)
{
  std::vector<GridTools::PeriodicFacePair<
    typename Triangulation<dim, spacedim>::cell_iterator>>
    periodicity_vector;

  GridTools::collect_periodic_faces(triangulation,
                                    bc_param.periodic_boundary_0,
                                    bc_param.periodic_boundary_1,
                                    bc_param.periodic_direction,
                                    periodicity_vector);
  triangulation.add_periodicity(periodicity_vector);
}

template void
read_mesh<2, 2>(const Parameters::Mesh                       &mesh_parameters,
                const bool                                    restart,
                const ConditionalOStream                     &pcout,
                parallel::DistributedTriangulationBase<2, 2> &triangulation,
                const Parameters::Lagrangian::BCDEM          &bc_params);

template void
read_mesh<2, 3>(const Parameters::Mesh                       &mesh_parameters,
                const bool                                    restart,
                const ConditionalOStream                     &pcout,
                parallel::DistributedTriangulationBase<2, 3> &triangulation,
                const Parameters::Lagrangian::BCDEM          &bc_params);

template void
read_mesh<3, 3>(const Parameters::Mesh                       &mesh_parameters,
                const bool                                    restart,
                const ConditionalOStream                     &pcout,
                parallel::DistributedTriangulationBase<3, 3> &triangulation,
                const Parameters::Lagrangian::BCDEM          &bc_params);
