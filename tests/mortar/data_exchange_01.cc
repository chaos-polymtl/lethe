// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: Create points on intersected mesh and determine owners on both
 * sides.
 */

#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/particle_handler.h>

#include <fstream>

// Lethe
#include <core/mortar_coupling_manager.h>

using namespace dealii;

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  const MPI_Comm comm = MPI_COMM_WORLD;

  const unsigned int dim                  = 2;
  const unsigned int n_global_refinements = 2;
  const unsigned int n_quadrature_points  = 3;
  const double       radius               = 1.0;
  const double       rotate               = 3.0;
  const double       rotate_pi            = 2 * numbers::PI * rotate / 360.0;

  // initialize triangulation
  parallel::distributed::Triangulation<dim> tria(comm);
  Triangulation<dim>                        tria_0, tria_1;

  // generate inner grid
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

  // generate outer grid
  GridGenerator::hyper_cube_with_cylindrical_hole(tria_1, radius, 2.0, true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  // create unique triangulation
  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 0, true, true);
  // store manifolds in merged triangulation
  tria.set_manifold(0, tria_0.get_manifold(0));
  tria.set_manifold(1, tria_0.get_manifold(1));
  tria.set_manifold(2, tria_1.get_manifold(0));

  tria.refine_global(n_global_refinements);

  const MortarManagerCircle<dim> mm(4 *
                                      Utilities::pow(2,
                                                     n_global_refinements + 1),
                                    radius,
                                    QGauss<dim>(n_quadrature_points),
                                    rotate_pi);

  const unsigned int n_points = mm.get_n_total_points();

  // convert local/ghost points to indices
  std::vector<double>                  local_values;
  std::vector<types::global_dof_index> is_local;
  std::vector<types::global_dof_index> is_ghost;

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == 0) || (face->boundary_id() == 2))
          {
            const auto indices = mm.get_indices(face->center());
            const auto points  = mm.get_points(face->center());

            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                unsigned int i = indices[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == 0)
                  {
                    id_local = i;
                    id_ghost = i + n_points;
                  }
                else if (face->boundary_id() == 2)
                  {
                    id_local = i + n_points;
                    id_ghost = i;
                  }

                is_local.emplace_back(id_local);
                is_ghost.emplace_back(id_ghost);

                for (unsigned int d = 0; d < dim; ++d)
                  local_values.emplace_back(points[ii][d]);
              }
          }

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  partitioner.reinit(is_local, is_ghost, comm);

  std::vector<double> ghost_values(local_values.size());
  partitioner.template export_to_ghosted_array<double, dim>(local_values,
                                                            ghost_values);

  for (unsigned int i = 0; i < local_values.size(); ++i)
    AssertThrow(std::abs(local_values[i] - ghost_values[i]) < 1e-8,
                ExcInternalError());

  ghost_values.assign(local_values.size(), 0.0);

  partitioner.template export_to_ghosted_array<double, 0>(local_values,
                                                          ghost_values,
                                                          dim);

  for (unsigned int i = 0; i < local_values.size(); ++i)
    AssertThrow(std::abs(local_values[i] - ghost_values[i]) < 1e-8,
                ExcInternalError());
}
