// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
}

template <int dim>
Point<dim>
rad_to_point(const double radius, const double rad)
{
  Point<dim> point;

  point[0] = radius * std::cos(rad);
  point[1] = radius * std::sin(rad);

  return point;
}

template <int dim>
class MortarManager
{
public:
  MortarManager(const unsigned int n_global_refinements,
                const unsigned int n_quadrature_points,
                const double       radius,
                const double       rotate)
    : n_global_refinements(n_global_refinements)
    , n_quadrature_points(n_quadrature_points)
    , radius(radius)
    , rotate(rotate)
    , quadrature(n_quadrature_points)
  {}

  std::pair<unsigned int, unsigned int>
  get_config(const double &rad) const
  {
    const double tolerance = 1e-8;
    const double delta =
      2 * numbers::PI / (4 * Utilities::pow(2, n_global_refinements + 1));
    const double rotate_pi = 2 * numbers::PI * rotate / 360.0;

    const double segment     = (rad - delta / 2) / delta;
    const double segment_rot = (rad - delta / 2 - rotate_pi) / delta;

    if (std::abs(rotate_pi / delta - std::round(rotate_pi / delta)) < tolerance)
      {
        // case 1: mesh is aligned
        return {0, std::round(segment)};
      }
    else
      {
        // case 2: mesh is not aligned
        if (std::abs(segment - std::round(segment)) < tolerance)
          return {2, std::round(segment)};
        else
          return {1, std::round(segment_rot)};
      }
  }

  unsigned int
  get_n_points() const
  {
    const auto [type, id] = get_config(0.0 /*not relevant*/);

    if (type == 0) // aligned
      {
        return 4 * Utilities::pow(2, n_global_refinements + 1) *
               n_quadrature_points;
      }
    else // inside/outside
      {
        return 8 * Utilities::pow(2, n_global_refinements + 1) *
               n_quadrature_points;
      }
  }

  std::vector<unsigned int>
  get_indices(const double &rad) const
  {
    const auto [type, id] = get_config(rad);

    if (type == 0) // aligned
      {
        std::vector<unsigned int> indices;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const unsigned int index = id * n_quadrature_points + q;

            AssertIndexRange(index, get_n_points());

            indices.emplace_back(index);
          }

        return indices;
      }
    else if (type == 1) // inside
      {
        std::vector<unsigned int> indices;

        for (unsigned int q = 0; q < n_quadrature_points * 2; ++q)
          {
            const unsigned int index =
              (id * n_quadrature_points * 2 + n_quadrature_points + q) %
              get_n_points();

            AssertIndexRange(index, get_n_points());

            indices.emplace_back(index);
          }

        return indices;
      }
    else // outside
      {
        std::vector<unsigned int> indices;

        for (unsigned int q = 0; q < n_quadrature_points * 2; ++q)
          {
            const unsigned int index = id * n_quadrature_points * 2 + q;

            AssertIndexRange(index, get_n_points());

            indices.emplace_back(index);
          }

        return indices;
      }
  }

  std::vector<Point<dim>>
  get_points(const double &rad) const
  {
    const auto [type, id] = get_config(rad);

    const double delta =
      2 * numbers::PI / (4 * Utilities::pow(2, n_global_refinements + 1));
    const double rotate_pi = 2 * numbers::PI * rotate / 360.0;

    if (type == 0) // aligned
      {
        std::vector<Point<dim>> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(
            rad_to_point<dim>(radius, (id + quadrature.point(q)[0]) * delta));

        return points;
      }
    else
      {
        double rad_0, rad_1, rad_2;

        double rot_min = rotate_pi - std::floor(rotate_pi / delta) * delta;

        if (type == 2)
          {
            rad_0 = id * delta;
            rad_1 = id * delta + rot_min;
            rad_2 = (id + 1) * delta;
          }
        else
          {
            rad_0 = id * delta + rot_min;
            rad_1 = (id + 1) * delta;
            rad_2 = (id + 1) * delta + rot_min;
          }

        std::vector<Point<dim>> points;

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(rad_to_point<dim>(radius,
                                                rad_0 + quadrature.point(q)[0] *
                                                          (rad_1 - rad_0)));

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          points.emplace_back(rad_to_point<dim>(radius,
                                                rad_1 + quadrature.point(q)[0] *
                                                          (rad_2 - rad_1)));

        return points;
      }
  }

private:
  const unsigned int n_global_refinements;
  const unsigned int n_quadrature_points;
  const double       radius;
  const double       rotate;
  QGauss<1>          quadrature;
};

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

  // create meshes
  parallel::distributed::Triangulation<dim> tria(comm);
  Triangulation<dim>                        tria_0, tria_1;

  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

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

  const MortarManager<dim> mm(n_global_refinements,
                              n_quadrature_points,
                              radius,
                              rotate);

  const unsigned int n_points = mm.get_n_points();

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
            const auto indices = mm.get_indices(point_to_rad(face->center()));

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
            const auto indices = mm.get_indices(point_to_rad(face->center()));

            // get points
            const auto points = mm.get_points(point_to_rad(face->center()));

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
