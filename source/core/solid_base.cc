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
 * Author: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020-
 */


#include <core/solid_base.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/data_out.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <memory.h>

#include <fstream>



template <int dim, int spacedim>
SolidBase<dim, spacedim>::SolidBase(
  std::shared_ptr<Parameters::NitscheSolid<spacedim>> &             param,
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>> fluid_tria,
  const unsigned int degree_velocity)
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
  , velocity(&param->solid_velocity)
  , degree_velocity(degree_velocity)
{}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::initial_setup()
{
  // initial_setup is called if the simulation is not a restarted one
  // Set-up of Nitsche triangulation, then particles (order important)
  setup_triangulation(false);
  setup_particles();
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_triangulation(const bool restart)
{
  if (param->solid_mesh.type == Parameters::Mesh::Type::gmsh)
    {
      // Grid creation
      GridIn<dim, spacedim> grid_in;
      // Attach triangulation to the grid
      grid_in.attach_triangulation(*solid_tria);
      // Read input gmsh file
      std::ifstream input_file(param->solid_mesh.file_name);
      grid_in.read_msh(input_file);
    }
  else if (param->solid_mesh.type == Parameters::Mesh::Type::dealii)
    {
      // deal.ii creates grid and attaches the solid triangulation
      GridGenerator::generate_from_name_and_arguments(
        *solid_tria,
        param->solid_mesh.grid_type,
        param->solid_mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - solid mesh will not be created");

  // Refine the solid triangulation to its initial size
  // NB: solid_tria should not be refined if loaded from a restart file
  // afterwards
  if (!restart)
    {
      solid_tria->refine_global(param->solid_mesh.initial_refinement);
      // Initialize dof handler for solid
      FE_Q<dim, spacedim> fe(1);
      solid_dh.distribute_dofs(fe);
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::load_triangulation(const std::string filename_tria)
{
  // Load solid triangulation from given file
  // TODO not functional for now, as not all information as passed with the load
  // function (see dealii documentation for classTriangulation) => change the
  // way the load works, or change the way the solid triangulation is handled
  std::ifstream in_folder(filename_tria.c_str());
  if (!in_folder)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename_tria + "> does not appear to exist!"));

  try
    {
      if (auto solid_tria =
            dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
              get_solid_triangulation().get()))
        {
          solid_tria->load(filename_tria.c_str());
        }
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "solid triangulation stored there."));
    }

  // Initialize dof handler for solid
  FE_Q<dim, spacedim> fe(1);
  solid_dh.distribute_dofs(fe);
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_particles_handler()
{
  unsigned int n_properties = (dim == 2 && spacedim == 3) ? 1 + spacedim : 1;

  solid_particle_handler =
    std::make_shared<Particles::ParticleHandler<spacedim>>();

  solid_particle_handler->initialize(*fluid_tria,
                                     StaticMappingQ1<spacedim>::mapping,
                                     n_properties);

  // Connect Nitsche particles to the fluid triangulation
  fluid_tria->signals.pre_distributed_refinement.connect(
    [&]() { solid_particle_handler->register_store_callback_function(); });
  fluid_tria->signals.post_distributed_refinement.connect(
    [&]() { solid_particle_handler->register_load_callback_function(false); });
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_particles()
{
  setup_particles_handler();

  // Initialize FEValues
  FE_Q<dim, spacedim> fe(1);

  QGauss<dim>                  quadrature(degree_velocity + 1);
  std::vector<Point<spacedim>> quadrature_points_vec;
  quadrature_points_vec.reserve(quadrature.size() *
                                solid_tria->n_locally_owned_active_cells());
  std::vector<std::vector<double>> properties;
  properties.reserve(quadrature.size() *
                     solid_tria->n_locally_owned_active_cells());

  UpdateFlags update_flags = update_JxW_values | update_quadrature_points;
  if (dim == 2 && spacedim == 3)
    update_flags = update_flags | update_normal_vectors;

  FEValues<dim, spacedim> fe_v(fe, quadrature, update_flags);

  // Fill quadrature points vector and properties, used to fill
  // solid_particle_handler
  for (const auto &cell : solid_dh.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_v.reinit(cell);
        const auto &points = fe_v.get_quadrature_points();
        const auto &JxW    = fe_v.get_JxW_values();
        if (dim == 2 && spacedim == 3)
          {
            const auto &normal_vectors = fe_v.get_normal_vectors();
            for (unsigned int q = 0; q < points.size(); ++q)
              {
                quadrature_points_vec.emplace_back(points[q]);
                std::vector<double> prop_i = {JxW[q],
                                              normal_vectors[q][0],
                                              normal_vectors[q][1]};
                if (spacedim == 3)
                  prop_i.push_back(normal_vectors[q][2]);
                properties.emplace_back(prop_i);
              }
          }
        else
          {
            for (unsigned int q = 0; q < points.size(); ++q)
              {
                quadrature_points_vec.emplace_back(points[q]);
                std::vector<double> prop_i = {JxW[q]};
                properties.emplace_back(prop_i);
              }
          }
      }

  // Compute fluid bounding box
  std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;

  // if Triangulation is a parallel::distributed::triangulation, use the naive
  // bounding box algorithm of deal.II
  if (auto tria =
        dynamic_cast<parallel::distributed::Triangulation<spacedim> *>(
          fluid_tria.get()))
    {
      // Increase number of bounding box level to ensure that insertion is
      // faster. Right now this is set to the maximum refinement level-1 which
      // is a decent heuristic.
      const unsigned int bounding_box_level = fluid_tria->n_global_levels() - 1;
      const auto         my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          *tria, IteratorFilters::LocallyOwnedCell(), bounding_box_level, true);
      global_fluid_bounding_boxes =
        Utilities::MPI::all_gather(mpi_communicator, my_bounding_box);
    }
  // else, use the more general boost rtree bounding boxes
  else
    {
      std::vector<BoundingBox<spacedim>> all_boxes;
      all_boxes.reserve(fluid_tria->n_locally_owned_active_cells());
      for (const auto cell : fluid_tria->active_cell_iterators())
        if (cell->is_locally_owned())
          all_boxes.emplace_back(cell->bounding_box());
      const auto tree        = pack_rtree(all_boxes);
      const auto local_boxes = extract_rtree_level(tree, 1);

      global_fluid_bounding_boxes =
        Utilities::MPI::all_gather(mpi_communicator, local_boxes);
    }

  // Fill solid particle handler
  solid_particle_handler->insert_global_particles(quadrature_points_vec,
                                                  global_fluid_bounding_boxes,
                                                  properties);

  // Number of particles used to assess that no particle is lost
  initial_number_of_particles = solid_particle_handler->n_global_particles();
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::load_particles(const std::string filename_part)
{
  // Setup particles handler
  setup_particles_handler();

  // Gather particle serialization information
  std::ifstream input(filename_part.c_str());
  AssertThrow(input, ExcFileNotOpen(filename_part));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> *solid_particle_handler;

  initial_number_of_particles = solid_particle_handler->n_global_particles();
}

template <int dim, int spacedim>
DoFHandler<dim, spacedim> &
SolidBase<dim, spacedim>::get_solid_dof_handler()
{
  return solid_dh;
}

template <int dim, int spacedim>
std::shared_ptr<Particles::ParticleHandler<spacedim>> &
SolidBase<dim, spacedim>::get_solid_particle_handler()
{
  return solid_particle_handler;
}

template <int dim, int spacedim>
std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>>
SolidBase<dim, spacedim>::get_solid_triangulation()
{
  return solid_tria;
}

template <int dim, int spacedim>
Function<spacedim> *
SolidBase<dim, spacedim>::get_solid_velocity()
{
  return velocity;
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::integrate_velocity(double time_step)
{
  const unsigned int sub_particles_iterations = param->particles_sub_iterations;
  AssertThrow(sub_particles_iterations >= 1,
              ExcMessage("Sub particles iterations must be 1 or larger"));
  double sub_iteration_relaxation = 1. / sub_particles_iterations;
  time_step                       = time_step * sub_iteration_relaxation;
  // Particle sub iterations divide the time step in a number of "sub
  // iterations". This allows the solver to use a larger CFL without
  // necessitating the use of the more complex particle location detection
  // approaches. The number of sub particles iterations must be chosen so that
  // the time_step / sub_particles_iterations  * velocity / cell size is smaller
  // than unity
  for (unsigned int it = 0; it < sub_particles_iterations; ++it)
    {
      for (auto particle = solid_particle_handler->begin();
           particle != solid_particle_handler->end();
           ++particle)
        {
          Point<spacedim> particle_location = particle->get_location();

          Tensor<1, spacedim> k1;
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k1[comp_i] = velocity->value(particle_location, comp_i);

          Point<spacedim>     p1 = particle_location + time_step / 2 * k1;
          Tensor<1, spacedim> k2;
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k2[comp_i] = velocity->value(p1, comp_i);

          Point<spacedim>     p2 = particle_location + time_step / 2 * k2;
          Tensor<1, spacedim> k3;
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k3[comp_i] = velocity->value(p2, comp_i);

          Point<spacedim>     p3 = particle_location + time_step * k3;
          Tensor<1, spacedim> k4;
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k4[comp_i] = velocity->value(p3, comp_i);

          particle_location += time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
          particle->set_location(particle_location);
        }
      solid_particle_handler->sort_particles_into_subdomains_and_cells();
    }

  if (initial_number_of_particles !=
      solid_particle_handler->n_global_particles())
    {
      if (this_mpi_process == 0)
        {
          std::cout << "Warning - Nitsche Particles have been lost"
                    << std::endl;
          std::cout << "Initial number of particles : "
                    << initial_number_of_particles << std::endl;
          std::cout << "Current number of particles : "
                    << solid_particle_handler->n_global_particles()
                    << std::endl;
        }
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::move_solid_triangulation(double time_step)
{
  const unsigned int n_dofs = solid_dh.n_dofs();
  std::vector<bool>  displacement(n_dofs, false);

  for (const auto &cell : solid_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0;
               i < GeometryInfo<spacedim>::vertices_per_cell;
               ++i)
            {
              if (!displacement[cell->vertex_index(i)])
                {
                  Point<spacedim> &vertex_position = cell->vertex(i);

                  Tensor<1, spacedim> k1;
                  for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                    k1[comp_i] = velocity->value(vertex_position, comp_i);

                  Point<spacedim>     p1 = vertex_position + time_step / 2 * k1;
                  Tensor<1, spacedim> k2;
                  for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                    k2[comp_i] = velocity->value(p1, comp_i);

                  Point<spacedim>     p2 = vertex_position + time_step / 2 * k2;
                  Tensor<1, spacedim> k3;
                  for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                    k3[comp_i] = velocity->value(p2, comp_i);

                  Point<spacedim>     p3 = vertex_position + time_step * k3;
                  Tensor<1, spacedim> k4;
                  for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                    k4[comp_i] = velocity->value(p3, comp_i);

                  vertex_position +=
                    time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

                  displacement[cell->vertex_index(i)] = true;
                }
            }
        }
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::print_particle_positions()
{
  for (auto particle = solid_particle_handler->begin();
       particle != solid_particle_handler->end();
       ++particle)
    {
      std::cout << "Particle " << particle->get_id() << " : "
                << particle->get_location() << std::endl;
    }
}

// Pre-compile the 2D, 3D and the 2D in 3D versions with the types that can
// occur
template class SolidBase<2>;
template class SolidBase<3>;
template class SolidBase<2, 3>;
