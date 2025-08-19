// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief MortarManagerCircle: output points as particles.
 */

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

#include "./tests.h"

using namespace dealii;

template <int dim, int spacedim>
void
output_mesh(const Triangulation<dim, spacedim> &tria,
            const unsigned int                  mapping_degree,
            const std::string                   file_name)
{
  MappingQ<dim, spacedim> mapping(mapping_degree);

  DataOut<dim, spacedim> data_out;

  if (dim == spacedim)
    {
      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);
    }

  data_out.attach_triangulation(tria);

  Vector<double> vector(tria.n_active_cells());
  vector = Utilities::MPI::this_mpi_process(tria.get_mpi_communicator());
  data_out.add_data_vector(vector, "ranks");

  data_out.build_patches(
    mapping,
    mapping_degree + 1,
    DataOut<dim, spacedim>::CurvedCellRegion::curved_inner_cells);

  data_out.write_vtu_in_parallel(file_name, tria.get_mpi_communicator());
}

template <int dim, int spacedim>
void
output_mesh(const std::vector<Point<spacedim>>     &points,
            const std::vector<std::vector<double>> &properties,
            const std::string                       file_name)
{
  BoundingBox<spacedim> bb(points);

  Triangulation<spacedim> tria;
  GridGenerator::hyper_rectangle(tria,
                                 bb.get_boundary_points().first,
                                 bb.get_boundary_points().second);

  MappingQ1<spacedim> mapping;

  std::vector<std::vector<BoundingBox<spacedim>>> bbs(1);
  bbs[0].emplace_back(bb);

  Particles::ParticleHandler<spacedim> particle_handler(tria, mapping, 2);
  particle_handler.insert_global_particles(points, bbs, properties);

  Particles::DataOut<spacedim> data_out_particles;
  data_out_particles.build_patches(
    particle_handler,
    {"0", "1"},
    {DataComponentInterpretation::component_is_scalar,
     DataComponentInterpretation::component_is_scalar});
  std::ofstream file_particles(file_name);
  data_out_particles.write_vtu(file_particles);
}

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

  // generated merged grid
  parallel::distributed::Triangulation<dim> tria(comm);
  hyper_cube_with_cylindrical_hole(radius, 2.0, rotate, tria);

  tria.refine_global(n_global_refinements);
  output_mesh<dim, dim>(tria, 3, "outer.0.vtu");

  const MortarManagerCircle<dim> mm(4 *
                                      Utilities::pow(2,
                                                     n_global_refinements + 1),
                                    radius,
                                    QGauss<dim>(n_quadrature_points),
                                    rotate_pi);

  const unsigned int n_points = mm.get_n_total_points();

  // convert local/ghost points to indices
  IndexSet is_local(n_points * 2);
  IndexSet is_ghost(n_points * 2);
  IndexSet is_points(n_points);

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == 0) || (face->boundary_id() == 2))
          {
            // get indices
            const auto mortar_indices =
              mm.get_mortar_indices(face->center(), face->boundary_id() == 0);

            std::vector<types::global_dof_index> indices(mm.get_n_points());

            for (unsigned int i = 0; i < indices.size(); ++i)
              indices[i] = (i < n_quadrature_points ? mortar_indices[0] :
                                                      mortar_indices[1]) *
                             n_quadrature_points +
                           (i % n_quadrature_points);

            for (const auto i : indices)
              {
                if (face->boundary_id() == 0)
                  {
                    is_local.add_index(i);
                    is_ghost.add_index(i + n_points);
                    is_points.add_index(i);
                  }
                else if (face->boundary_id() == 2)
                  {
                    is_local.add_index(i + n_points);
                    is_ghost.add_index(i);
                    is_points.add_index(i);
                  }
              }
          }

  is_ghost.subtract_set(is_local);

  // determine owner of indices
  const auto ghost_owners =
    Utilities::MPI::compute_index_owner(is_local, is_ghost, comm);


  // output result
  std::vector<Point<dim>>          relevant_points(is_points.n_elements());
  std::vector<std::vector<double>> properties(is_points.n_elements(),
                                              std::vector<double>(2, -1));

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == 0) || (face->boundary_id() == 2))
          {
            // get indices
            const auto mortar_indices =
              mm.get_mortar_indices(face->center(), face->boundary_id() == 0);

            std::vector<types::global_dof_index> indices(mm.get_n_points());

            for (unsigned int i = 0; i < indices.size(); ++i)
              indices[i] = (i < n_quadrature_points ? mortar_indices[0] :
                                                      mortar_indices[1]) *
                             n_quadrature_points +
                           (i % n_quadrature_points);

            // get points
            const auto points =
              mm.get_points(face->center(), face->boundary_id() == 0);

            for (unsigned int i = 0; i < indices.size(); ++i)
              {
                const auto index =
                  is_points.index_within_set(indices[i]) % n_points;

                relevant_points[index] = points[i];

                if (face->boundary_id() == 0)
                  properties[index][0] =
                    is_local.is_element(indices[i] + n_points) ?
                      Utilities::MPI::this_mpi_process(comm) :
                      ghost_owners[is_ghost.index_within_set(indices[i] +
                                                             n_points)];
                else if (face->boundary_id() == 2)
                  properties[index][1] =
                    is_local.is_element(indices[i]) ?
                      Utilities::MPI::this_mpi_process(comm) :
                      ghost_owners[is_ghost.index_within_set(indices[i])];
              }
          }

  output_mesh<dim, dim>(
    relevant_points,
    properties,
    "points" + std::to_string(Utilities::MPI::this_mpi_process(comm)) + ".vtu");
}
