/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2019-
 */
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/particles/data_out.h>

#include <core/grids.h>
#include <core/parameters.h>
#include <solvers/solid_base.h>


template <int dim, int spacedim>
SolidBase<dim, spacedim>::SolidBase(
  Parameters::Nitsche                                                &param,
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>>   fluid_tria,
  const unsigned int                                                  degree_velocity)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , solid_tria(
      dynamic_cast<parallel::DistributedTriangulationBase<dim, spacedim> *>(
        new parallel::distributed::Triangulation<dim, spacedim>(
          mpi_communicator,
          typename Triangulation<dim, spacedim>::MeshSmoothing(
            Triangulation<dim, spacedim>::smoothing_on_refinement |
            Triangulation<dim, spacedim>::smoothing_on_coarsening))))
  , fluid_tria(fluid_tria)
  , solid_dh(*solid_tria)
  , param(param)
  , degree_velocity(degree_velocity)
{}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::initial_setup()
{
  FE_Nothing<dim, spacedim> solid_fe;
  solid_dh.distribute_dofs(solid_fe);

  if (param.solid_mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim, spacedim> grid_in;
      grid_in.attach_triangulation(*solid_tria);
      std::ifstream input_file(param.solid_mesh.file_name);
      grid_in.read_msh(input_file);
    }
  else if (param.solid_mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        *solid_tria,
        param.solid_mesh.grid_type,
        param.solid_mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - solid mesh will not be created");
}


template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_particles()
{
  QGauss<dim>        quadrature(degree_velocity + 1);
  const unsigned int n_properties = 1;
  solid_particle_handler->initialize(*fluid_tria,
                                     StaticMappingQ1<spacedim>::mapping,
                                     n_properties);
  std::vector<Point<spacedim>> quadrature_points_vec;
  quadrature_points_vec.reserve(quadrature.size() *
                                solid_tria->n_locally_owned_active_cells());
  std::vector<std::vector<double>> properties;
  properties.reserve(quadrature.size() *
                     solid_tria->n_locally_owned_active_cells());

  FE_Q<dim, spacedim>     fe(1);
  FEValues<dim, spacedim> fe_v(fe,
                               quadrature,
                               update_JxW_values | update_quadrature_points);
  for (const auto &cell : solid_dh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_v.reinit(cell);
        const auto &points = fe_v.get_quadrature_points();
        const auto &JxW    = fe_v.get_JxW_values();
        for (unsigned int q = 0; q < points.size(); ++q)
          {
            quadrature_points_vec.emplace_back(points[q]);
            properties.emplace_back(std::vector<double>(n_properties, JxW[q]));
          }
      }

  std::vector<BoundingBox<spacedim>> all_boxes;
  all_boxes.reserve(fluid_tria->n_locally_owned_active_cells());
  for (const auto cell : fluid_tria->active_cell_iterators())
    if (cell->is_locally_owned())
      all_boxes.emplace_back(cell->bounding_box());
  const auto tree        = pack_rtree(all_boxes);
  const auto local_boxes = extract_rtree_level(tree, 1);

  std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;
  global_fluid_bounding_boxes =
    Utilities::MPI::all_gather(mpi_communicator, local_boxes);

  solid_particle_handler->insert_global_particles(quadrature_points_vec,
                                                  global_fluid_bounding_boxes,
                                                  properties);

  fluid_tria->signals.pre_distributed_refinement.connect(
    [&]() { solid_particle_handler->register_store_callback_function(); });
  fluid_tria->signals.post_distributed_refinement.connect(
    [&]() { solid_particle_handler->register_load_callback_function(false); });

    setup_done = true;
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::output_particles(std::string fprefix) const
{
  Particles::DataOut<spacedim, spacedim> particles_out;
  particles_out.build_patches(*solid_particle_handler);
  const std::string filename = (fprefix + ".vtu");
  particles_out.write_vtu_in_parallel("/" + filename, mpi_communicator);
}

template <int dim, int spacedim>
std::shared_ptr<Particles::ParticleHandler<spacedim>>
SolidBase<dim, spacedim>::get_solid_particle_handler()
{
  if (!setup_done) {
    initial_setup();
    setup_particles();
  }
  return solid_particle_handler;
}

// Pre-compile the 2D, 3D and the 2D in 3D versions with the types that can
// occur
template class SolidBase<2>;
template class SolidBase<3>;
template class SolidBase<2, 3>;
