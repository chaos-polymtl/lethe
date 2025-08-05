// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: Collect indices and create sparsity pattern.
 */

#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
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
  const unsigned int fe_degree            = 1;
  const unsigned int n_quadrature_points  = 1;
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

  FE_Q<dim> fe(fe_degree);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const MortarManagerCircle<dim> mm(4 *
                                      Utilities::pow(2,
                                                     n_global_refinements + 1),
                                    radius,
                                    QGauss<dim>(n_quadrature_points),
                                    rotate_pi);

  const unsigned int n_points = mm.get_n_total_points();

  // convert local/ghost points to indices
  std::vector<types::global_dof_index> local_values;
  std::vector<types::global_dof_index> is_local;
  std::vector<types::global_dof_index> is_ghost;

  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == 0) || (face->boundary_id() == 2))
          {
            const auto indices = mm.get_indices(face->center());

            std::vector<types::global_dof_index> local_dofs(
              fe.n_dofs_per_cell());
            cell->get_dof_indices(local_dofs);

            for (const auto i : indices)
              {
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

                for (const auto i : local_dofs)
                  local_values.emplace_back(i);
              }
          }

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  partitioner.reinit(is_local, is_ghost, comm);

  std::vector<types::global_dof_index> ghost_values(local_values.size());

  partitioner.template export_to_ghosted_array<types::global_dof_index, 0>(
    local_values, ghost_values, fe.n_dofs_per_cell());

  if (Utilities::MPI::n_mpi_processes(comm) == 1)
    {
      // create sparsity pattern
      const unsigned int n_dofs_per_cell = fe.n_dofs_per_cell();

      DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                      dof_handler.n_dofs());

      // 1) cell entries
      std::vector<types::global_dof_index> dof_indices;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          dof_indices.resize(cell->get_fe().n_dofs_per_cell());
          cell->get_dof_indices(dof_indices);

          for (const auto i : dof_indices)
            dynamic_sparsity_pattern.add_row_entries(i, dof_indices);
        }

      // 2) coupling entries
      for (unsigned int i = 0; i < local_values.size(); i += n_dofs_per_cell)
        {
          std::vector<types::global_dof_index> a(local_values.begin() + i,
                                                 local_values.begin() + i +
                                                   n_dofs_per_cell);
          std::vector<types::global_dof_index> b(ghost_values.begin() + i,
                                                 ghost_values.begin() + i +
                                                   n_dofs_per_cell);

          for (const auto i : a)
            dynamic_sparsity_pattern.add_row_entries(i, b);
          for (const auto i : b)
            dynamic_sparsity_pattern.add_row_entries(i, a);
        }

      // 3) write sparsity pattern to file
      SparsityPattern sparsity_pattern;
      sparsity_pattern.copy_from(dynamic_sparsity_pattern);
      std::ofstream out("temp.svg");
      sparsity_pattern.print_svg(out);
    }
}
