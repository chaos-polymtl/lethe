/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2024 by the Lethe authors
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
 * Mortar: Create points on intersected mesh and determine owners on both sides.
 *
 * ---------------------------------------------------------------------*/

#include <deal.II/base/floating_point_comparator.h>

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
  vector = Utilities::MPI::this_mpi_process(tria.get_communicator());
  data_out.add_data_vector(vector, "ranks");

  data_out.build_patches(
    mapping,
    mapping_degree + 1,
    DataOut<dim, spacedim>::CurvedCellRegion::curved_inner_cells);

  data_out.write_vtu_in_parallel(file_name, tria.get_communicator());
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



template <int dim, int spacedim>
std::tuple<std::vector<Point<spacedim>>, std::vector<double>>
compute_quadrature(double             rad_00,
                   double             rad_01,
                   double             rad_10,
                   double             rad_11,
                   const double       radius,
                   const unsigned int n_quadrature_points,
                   const unsigned int mapping_degree = 10)
{
  if (rad_00 > rad_10) // normalize
    {
      std::swap(rad_00, rad_10);
      std::swap(rad_01, rad_11);
    }

  if (rad_01 <= rad_10)
    return {}; // no cut

  const double left_rad  = std::max(rad_00, rad_10);
  const double right_rad = std::min(rad_01, rad_11);

  std::vector<Point<spacedim>> points;
  points.emplace_back(std::cos(left_rad) * radius, std::sin(left_rad) * radius);
  points.emplace_back(std::cos(right_rad) * radius,
                      std::sin(right_rad) * radius);

  std::vector<CellData<dim>> cells;
  CellData<dim>              cell(2);
  cell.vertices[0] = 0;
  cell.vertices[1] = 1;
  cells.emplace_back(cell);

  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(points, cells, {});

  tria.set_manifold(0, SphericalManifold<dim, spacedim>());
  tria.set_all_manifold_ids(0);

  MappingQ<dim, spacedim>   mapping(mapping_degree);
  QGauss<dim>               quadrature(n_quadrature_points);
  FE_Nothing<dim, spacedim> fe;

  FEValues<dim, spacedim> phi(mapping,
                              fe,
                              quadrature,
                              update_quadrature_points | update_JxW_values);
  phi.reinit(tria.begin());

  return {phi.get_quadrature_points(), phi.get_JxW_values()};
}


template <int dim>
double
point_to_rad(const Point<dim> &point)
{
  const double temp = std::atan(std::abs(point[1]) / std::abs(point[0]));

  if (point[1] >= 0.0)
    {
      if (point[0] >= 0.0)
        return temp;
      else
        return numbers::PI - temp;
    }
  else
    {
      if (point[0] >= 0.0)
        return 2 * numbers::PI - temp;
      else
        return numbers::PI + temp;
    }
};

template <int structdim, int dim, int spacedim>
std::vector<Point<spacedim>>
compute_quadrature(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &cell_0,
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &cell_1,
  const double                                                radius,
  const unsigned int n_quadrature_points,
  const unsigned int mapping_degree = 10)
{
  AssertDimension(structdim, 1);
  AssertDimension(spacedim, 2);

  // helper function to 1) guarantee that all segments have the same orientation
  // and 2) periodicities are handled correctly -> split up into two segment
  const auto create_sections = [](const auto &cell) {
    const auto cross_product = [](const auto &p0, const auto &p1) {
      Tensor<1, 3, double> t0;
      Tensor<1, 3, double> t1;

      for (unsigned int i = 0; i < 2; ++i)
        {
          t0[i] = p0[i];
          t1[i] = p1[i];
        }

      return cross_product_3d(t0, t1)[2];
    };

    double rad_0 = point_to_rad(cell->vertex(0));
    double rad_1 = point_to_rad(cell->vertex(1));

    if (cross_product(cell->vertex(0), cell->vertex(1)) < 0.0)
      std::swap(rad_0, rad_1);

    std::vector<std::pair<double, double>> sections;

    if (rad_0 < rad_1) // normal
      {
        sections.emplace_back(rad_0, rad_1);
      }
    else // periodic
      {
        sections.emplace_back(rad_0, 2 * numbers::PI);
        sections.emplace_back(0, rad_1);
      }

    return sections;
  };

  const auto subsections_0 = create_sections(cell_0);
  const auto subsections_1 = create_sections(cell_1);

  std::vector<Point<spacedim>> all_points;

  for (const auto &subsection_0 : subsections_0)
    for (const auto &subsection_1 : subsections_1)
      {
        const auto [points, JxWs] =
          compute_quadrature<structdim, spacedim>(subsection_0.first,
                                                  subsection_0.second,
                                                  subsection_1.first,
                                                  subsection_1.second,
                                                  radius,
                                                  n_quadrature_points,
                                                  mapping_degree);

        all_points.insert(all_points.end(), points.begin(), points.end());
      }

  return all_points;
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

  parallel::distributed::Triangulation<dim> tria(comm);
  Triangulation<dim>                        tria_0, tria_1;
  Triangulation<dim - 1, dim>               tria_0_surface, tria_1_surface;

  // create meshes
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(3, tria_0);

  GridGenerator::hyper_cube_with_cylindrical_hole(tria_1, radius, 2.0, true);
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 0, true, true);

  tria.set_manifold(0, tria_0.get_manifold(0));
  tria.set_manifold(1, tria_0.get_manifold(1));
  tria.set_manifold(2, tria_1.get_manifold(0));

  tria.refine_global(n_global_refinements);
  output_mesh<dim, dim>(tria, 3, "outer.0.vtu");

  // create surface meshes
  GridGenerator::hyper_sphere(tria_0_surface, {}, radius);
  GridTools::rotate(3, tria_0_surface);
  tria_0_surface.refine_global(n_global_refinements + 1);

  GridGenerator::hyper_sphere(tria_1_surface, {}, radius);
  tria_1_surface.refine_global(n_global_refinements + 1);

  // collect all surface points -> oracle
  std::vector<Point<dim>> all_points;

  const FloatingPointComparator<double> comparator(1e-6);

  std::map<double, std::vector<unsigned int>, FloatingPointComparator<double>>
    all_points_0(comparator);
  std::map<double, std::vector<unsigned int>, FloatingPointComparator<double>>
    all_points_1(comparator);

  for (const auto &face_0 : tria_0_surface.active_cell_iterators())
    for (const auto &face_1 : tria_1_surface.active_cell_iterators())
      {
        const auto points = compute_quadrature<dim - 1, dim - 1, dim>(
          face_0, face_1, radius, n_quadrature_points);

        for (const auto &p : points)
          {
            all_points_0[point_to_rad(face_0->center())].emplace_back(
              all_points.size());
            all_points_1[point_to_rad(face_1->center())].emplace_back(
              all_points.size());
            all_points.emplace_back(p);
          }
      }

  // convert local/ghost points to indices
  IndexSet is_local(all_points.size() * 2);
  IndexSet is_ghost(all_points.size() * 2);
  IndexSet is_points(all_points.size());

  for (const auto &cell_0 : tria.active_cell_iterators())
    if (cell_0->is_locally_owned())
      for (const auto &face_0 : cell_0->face_iterators())
        if (face_0->boundary_id() == 0)
          {
            const auto &indices = all_points_0[point_to_rad(face_0->center())];

            for (const auto i : indices)
              {
                is_local.add_index(i + 0);
                is_ghost.add_index(i + all_points.size());
                is_points.add_index(i + 0);
              }
          }

  for (const auto &cell_1 : tria.active_cell_iterators())
    if (cell_1->is_locally_owned())
      for (const auto &face_1 : cell_1->face_iterators())
        if (face_1->boundary_id() == 2)
          {
            const auto &indices = all_points_1[point_to_rad(face_1->center())];

            for (const auto i : indices)
              {
                is_local.add_index(i + all_points.size());
                is_ghost.add_index(i + 0);
                is_points.add_index(i + 0);
              }
          }

  is_ghost.subtract_set(is_local);

  // determine owner of indices
  const auto ghost_owners =
    Utilities::MPI::compute_index_owner(is_local, is_ghost, comm);


  // output result
  std::vector<Point<dim>>          relevant_points(is_points.n_elements());
  std::vector<std::vector<double>> properties(is_points.n_elements(),
                                              std::vector<double>(2));

  for (const auto &cell : tria.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == 0) || (face->boundary_id() == 2))
          {
            // get indices
            const auto &indices = (face->boundary_id() == 0) ?
                                    all_points_0[point_to_rad(face->center())] :
                                    all_points_1[point_to_rad(face->center())];

            // get points
            std::vector<Point<dim>> points(indices.size());
            for (unsigned int i = 0; i < indices.size(); ++i)
              points[i] = all_points[indices[i]];

            for (unsigned int i = 0; i < indices.size(); ++i)
              {
                const auto index =
                  is_points.index_within_set(indices[i]) % all_points.size();

                relevant_points[index] = points[i];

                if (face->boundary_id() == 0)
                  properties[index][0] =
                    is_local.is_element(indices[i] + all_points.size()) ?
                      Utilities::MPI::this_mpi_process(comm) :
                      ghost_owners[is_ghost.index_within_set(
                        indices[i] + all_points.size())];
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
