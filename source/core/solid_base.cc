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

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>



template <int dim, int spacedim>
SolidBase<dim, spacedim>::SolidBase(
  std::shared_ptr<Parameters::NitscheObject<spacedim>> &            param,
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>> fluid_tria,
  std::shared_ptr<Mapping<spacedim>> fluid_mapping)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , fluid_tria(fluid_tria)
  , fluid_mapping(fluid_mapping)
  , param(param)
  , velocity(&param->solid_velocity)
  , temperature(&param->solid_temperature)
{
  if (param->solid_mesh.simplex)
    {
      // for simplex meshes
      fe            = std::make_shared<FE_SimplexP<dim, spacedim>>(1);
      solid_mapping = std::make_shared<MappingFE<dim, spacedim>>(*fe);
      quadrature =
        std::make_shared<QGaussSimplex<dim>>(param->number_quadrature_points);
      solid_tria = std::make_shared<
        parallel::fullydistributed::Triangulation<dim, spacedim>>(
        this->mpi_communicator);

      displacement_fe =
        std::make_shared<FESystem<dim, spacedim>>(FE_SimplexP<dim, spacedim>(1),
                                                  spacedim);
    }
  else
    {
      // Usual case, for quad/hex meshes
      fe            = std::make_shared<FE_Q<dim, spacedim>>(1);
      solid_mapping = std::make_shared<MappingQGeneric<dim, spacedim>>(1);
      quadrature =
        std::make_shared<QGauss<dim>>(param->number_quadrature_points);
      solid_tria =
        std::make_shared<parallel::distributed::Triangulation<dim, spacedim>>(
          mpi_communicator,
          typename Triangulation<dim, spacedim>::MeshSmoothing(
            Triangulation<dim, spacedim>::smoothing_on_refinement |
            Triangulation<dim, spacedim>::smoothing_on_coarsening));

      displacement_fe =
        std::make_shared<FESystem<dim, spacedim>>(FE_Q<dim, spacedim>(1),
                                                  spacedim);
    }

  // Dof handler associated with particles
  solid_dh.clear();
  solid_dh.reinit(*this->solid_tria);

  // Dof handler associated with mesh displacement
  displacement_dh.clear();
  displacement_dh.reinit(*this->solid_tria);
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::initial_setup()
{
  // initial_setup is called if the simulation is not a restarted one
  // Set-up of Nitsche triangulation, then particles (order important)
  setup_triangulation(false);
  setup_displacement();
  setup_particles();
}



template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_triangulation(const bool restart)
{
  if (param->solid_mesh.type == Parameters::Mesh::Type::gmsh)
    {
      if (param->solid_mesh.simplex)
        {
          auto        comm      = solid_tria->get_communicator();
          std::string file_name = param->solid_mesh.file_name;

          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation_in_groups<dim, spacedim>(
              [file_name](dealii::Triangulation<dim, spacedim> &basetria) {
                GridIn<dim, spacedim> grid_in;
                grid_in.attach_triangulation(basetria);
                std::ifstream input_file(file_name);

                grid_in.read_msh(input_file);
              },
              [](dealii::Triangulation<dim, spacedim> &basetria,
                 const MPI_Comm                        comm,
                 const unsigned int /*group_size*/) {
                GridTools::partition_triangulation(
                  Utilities::MPI::n_mpi_processes(comm), basetria);
              },
              comm,
              Utilities::MPI::n_mpi_processes(comm) /* group size */,
              dealii::Triangulation<dim, spacedim>::none);

          solid_tria->create_triangulation(construction_data);
        }
      else
        { // Grid creation
          GridIn<dim, spacedim> grid_in;
          // Attach triangulation to the grid
          grid_in.attach_triangulation(*solid_tria);
          // Read input gmsh file
          std::ifstream input_file(param->solid_mesh.file_name);
          grid_in.read_msh(input_file);
        }
    }
  else if (param->solid_mesh.type == Parameters::Mesh::Type::dealii)
    {
      if (param->solid_mesh.simplex)
        { // TODO Using dealii generated meshes with simplices in Nitsche solver
          // generates and error. "boundary_id !=
          // numbers::internal_face_boundary_id"
          Triangulation<dim, spacedim> temporary_quad_triangulation;
          GridGenerator::generate_from_name_and_arguments(
            temporary_quad_triangulation,
            param->solid_mesh.grid_type,
            param->solid_mesh.grid_arguments);

          // initial refinement
          const int initial_refinement = param->solid_mesh.initial_refinement;
          temporary_quad_triangulation.refine_global(initial_refinement);
          // flatten the triangulation
          Triangulation<dim, spacedim> flat_temp_quad_triangulation;
          GridGenerator::flatten_triangulation(temporary_quad_triangulation,
                                               flat_temp_quad_triangulation);

          Triangulation<dim, spacedim> temporary_tri_triangulation(
            Triangulation<dim, spacedim>::limit_level_difference_at_vertices);
          GridGenerator::convert_hypercube_to_simplex_mesh(
            flat_temp_quad_triangulation, temporary_tri_triangulation);

          GridTools::partition_triangulation_zorder(
            Utilities::MPI::n_mpi_processes(solid_tria->get_communicator()),
            temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              solid_tria->get_communicator(),
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          solid_tria->create_triangulation(construction_data);
        }
      else
        { // deal.ii creates grid and attaches the solid triangulation
          GridGenerator::generate_from_name_and_arguments(
            *solid_tria,
            param->solid_mesh.grid_type,
            param->solid_mesh.grid_arguments);
        }
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - solid mesh will not be created");

  // Rotate the triangulation
  rotate_grid(param->solid_mesh.rotation_angle,
              param->solid_mesh.rotation_axis);

  // Translate the triangulation
  translate_grid(param->solid_mesh.translation);
  // Refine the solid triangulation to its initial size
  // NB: solid_tria should not be refined if loaded from a restart file
  // afterwards
  if (!restart)
    {
      if (param->solid_mesh.simplex)
        {
          // Simplex triangulation refinement isn't possible yet
        }
      else
        {
          solid_tria->refine_global(param->solid_mesh.initial_refinement);
        }
      solid_dh.distribute_dofs(*fe);
    }
}

template <>
void
SolidBase<2, 2>::rotate_grid(const double angle, const Tensor<1, 3> /*axis*/)
{
  GridTools::rotate(angle, *solid_tria);
}
template <>
void
SolidBase<2, 3>::rotate_grid(const double angle, const Tensor<1, 3> axis)
{
  GridTools::rotate(axis, angle, *solid_tria);
}
template <>
void
SolidBase<3, 3>::rotate_grid(const double angle, const Tensor<1, 3> axis)
{
  GridTools::rotate(axis, angle, *solid_tria);
}

template <>
void
SolidBase<2, 2>::translate_grid(const Tensor<1, 3> translation)
{
  GridTools::shift(Tensor<1, 2>({translation[0], translation[1]}), *solid_tria);
}

template <>
void
SolidBase<2, 3>::translate_grid(const Tensor<1, 3> translation)
{
  GridTools::shift(translation, *solid_tria);
}

template <>
void
SolidBase<3, 3>::translate_grid(const Tensor<1, 3> translation)
{
  GridTools::shift(translation, *solid_tria);
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
  solid_dh.distribute_dofs(*fe);
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_particles_handler()
{
  unsigned int n_properties = (dim == 2 && spacedim == 3) ? 1 + spacedim : 1;

  solid_particle_handler =
    std::make_shared<Particles::ParticleHandler<spacedim>>();

  // Put the proper triangulation and mapping for more general cases
  solid_particle_handler->initialize(*fluid_tria, *fluid_mapping, n_properties);
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_particles()
{
  setup_particles_handler();

  // Initialize FEValues
  std::vector<Point<spacedim>> quadrature_points_vec;
  quadrature_points_vec.reserve(quadrature->size() *
                                solid_tria->n_locally_owned_active_cells());
  unsigned int n_properties = (dim == 2 && spacedim == 3) ? 1 + spacedim : 1;

  std::vector<std::vector<double>> properties;
  properties.reserve(quadrature->size() *
                     solid_tria->n_locally_owned_active_cells() * n_properties);

  UpdateFlags update_flags = update_JxW_values | update_quadrature_points;
  if (dim == 2 && spacedim == 3)
    update_flags = update_flags | update_normal_vectors;

  FEValues<dim, spacedim> fe_v(*solid_mapping, *fe, *quadrature, update_flags);

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



  // Use the more general boost rtree bounding boxes
  std::vector<BoundingBox<spacedim>> all_boxes;
  all_boxes.reserve(fluid_tria->n_locally_owned_active_cells());
  for (const auto &cell : fluid_tria->active_cell_iterators())
    if (cell->is_locally_owned())
      all_boxes.emplace_back(cell->bounding_box());
  const auto tree        = pack_rtree(all_boxes);
  const auto local_boxes = extract_rtree_level(
    tree, std::max(int(log10(fluid_tria->n_locally_owned_active_cells())), 1));

  global_fluid_bounding_boxes =
    Utilities::MPI::all_gather(mpi_communicator, local_boxes);

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
DoFHandler<dim, spacedim> &
SolidBase<dim, spacedim>::get_displacement_dof_handler()
{
  return displacement_dh;
}


template <int dim, int spacedim>
TrilinosWrappers::MPI::Vector &
SolidBase<dim, spacedim>::get_displacement_vector()
{
  displacement_relevant = displacement;
  return displacement_relevant;
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
Function<spacedim> *
SolidBase<dim, spacedim>::get_solid_temperature()
{
  return temperature;
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::update_temperature_time(double time)
{
  temperature->set_time(time);
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::integrate_velocity(double time_step,
                                             double initial_time)
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
          // Calculate the next position of the particle using the RK4
          // time integration scheme
          // x(t+dt) = x(t) + 1/6 (k1+2*k2+2*k3+k4)
          // k1 = dt * v(t,x(t))
          // k2 = dt * v(t+0.5*dt,x+k1/2)
          // k3 = dt * v(t+0.5*dt,x+k2/2)
          // k4 = dt * vt(t+dt,x+k3)
          // The four stages (1 to 4) are built successively
          // For each stage, the time of the function must be "reset"
          // to ensure adequate evaluation of the RK4 scheme.
          // See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods for
          // more details

          Point<spacedim> particle_location = particle->get_location();

          Tensor<1, spacedim> k1;
          velocity->set_time(initial_time);
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k1[comp_i] = velocity->value(particle_location, comp_i);

          Point<spacedim>     p1 = particle_location + time_step / 2 * k1;
          Tensor<1, spacedim> k2;
          velocity->set_time(initial_time + time_step / 2);
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k2[comp_i] = velocity->value(p1, comp_i);

          Point<spacedim>     p2 = particle_location + time_step / 2 * k2;
          Tensor<1, spacedim> k3;
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k3[comp_i] = velocity->value(p2, comp_i);

          Point<spacedim>     p3 = particle_location + time_step * k3;
          Tensor<1, spacedim> k4;
          velocity->set_time(initial_time + time_step);
          for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
            k4[comp_i] = velocity->value(p3, comp_i);

          particle_location += time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
          particle->set_location(particle_location);
        }
      solid_particle_handler->sort_particles_into_subdomains_and_cells();
      initial_time += time_step;
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
          std::cout << "Try increasing the number of particles sub iterations"
                    << std::endl;

          // Stops the simulation if particles have been lost
          // Future improvement would be to dynamically increase the number of
          // particles subiteration if particles have been lost
          if (param->stop_particles_lost)
            {
              std::cerr
                << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
              std::cerr
                << "Exception on processing: " << std::endl
                << "Nitsche Particles have been lost" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
              std::exit(1);
            }
        }
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::move_solid_triangulation(const double time_step,
                                                   const double initial_time)
{
  // First we calculate the displacement that the solid triangulation
  // will feel during this iteration.
  // Then, we apply this displacement to the triangulation
  const unsigned int      dofs_per_cell = displacement_fe->dofs_per_cell;
  const std::vector<bool> locally_owned_vertices =
    GridTools::get_locally_owned_vertices(*this->solid_tria);
  std::vector<bool> vertex_moved(this->solid_tria->n_vertices(), false);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


  for (const auto &cell : displacement_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int vertex = 0;
               vertex < cell->reference_cell().n_vertices();
               ++vertex)
            {
              // Calculate the next position of the particle using the RK4
              // time integration scheme
              // x(t+dt) = x(t) + 1/6 (k1+2*k2+2*k3+k4)
              // k1 = dt * v(t,x(t))
              // k2 = dt * v(t+0.5*dt,x+k1/2)
              // k3 = dt * v(t+0.5*dt,x+k2/2)
              // k4 = dt * vt(t+dt,x+k3)
              // The four stages (1 to 4) are built successively
              // For each stage, the time of the function must be "reset"
              // to ensure adequate evaluation of the RK4 scheme.
              // See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
              // for more details

              const auto       dof_index = cell->vertex_dof_index(vertex, 0);
              Point<spacedim> &vertex_position  = cell->vertex(vertex);
              const unsigned   global_vertex_no = cell->vertex_index(vertex);
              if (vertex_moved[global_vertex_no] ||
                  !locally_owned_vertices[global_vertex_no])
                continue;

              Tensor<1, spacedim> k1;
              velocity->set_time(initial_time);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k1[comp_i] = velocity->value(vertex_position, comp_i);

              Point<spacedim>     p1 = vertex_position + time_step / 2 * k1;
              Tensor<1, spacedim> k2;
              velocity->set_time(initial_time + time_step / 2);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k2[comp_i] = velocity->value(p1, comp_i);

              Point<spacedim>     p2 = vertex_position + time_step / 2 * k2;
              Tensor<1, spacedim> k3;
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k3[comp_i] = velocity->value(p2, comp_i);

              Point<spacedim>     p3 = vertex_position + time_step * k3;
              Tensor<1, spacedim> k4;
              velocity->set_time(initial_time + time_step);
              for (unsigned int comp_i = 0; comp_i < spacedim; ++comp_i)
                k4[comp_i] = velocity->value(p3, comp_i);

              auto vertex_displacement =
                time_step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);

              vertex_position += vertex_displacement;

              for (unsigned d = 0; d < spacedim; ++d)
                displacement[dof_index + d] =
                  displacement[dof_index + d] + vertex_displacement[d];

              vertex_moved[global_vertex_no] = true;
            }
        }
    }
  this->solid_tria->communicate_locally_moved_vertices(locally_owned_vertices);
}


template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::move_solid_triangulation_with_displacement()
{
  const unsigned int     dofs_per_cell = displacement_fe->dofs_per_cell;
  std::set<unsigned int> dof_vertex_displaced;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : displacement_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int vertex = 0;
               vertex < cell->reference_cell().n_vertices();
               ++vertex)
            {
              if (!dof_vertex_displaced.count(cell->vertex_index(vertex)))
                {
                  const auto dof_index = cell->vertex_dof_index(vertex, 0);
                  Point<spacedim> &vertex_position = cell->vertex(vertex);

                  for (unsigned d = 0; d < spacedim; ++d)
                    {
                      vertex_position[d] = vertex_position[d] +
                                           displacement_relevant[dof_index + d];
                    }

                  dof_vertex_displaced.insert(cell->vertex_index(vertex));
                }
            }
        }
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::print_particle_positions()
{
  std::map<int, Particles::ParticleIterator<spacedim>> local_particles;
  for (auto particle = solid_particle_handler->begin();
       particle != solid_particle_handler->end();
       ++particle)
    {
      local_particles.insert({particle->get_id(), particle});
    }

  for (auto &iterator : local_particles)
    {
      unsigned int id                = iterator.first;
      auto         particle          = iterator.second;
      auto         particle_location = particle->get_location();

      std::cout << std::fixed << std::setprecision(0) << id << " "
                << std::setprecision(2) << particle_location << std::endl;
    }
}


template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::write_checkpoint(std::string prefix)
{
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::Vector,
                                          DoFHandler<dim, spacedim>>
    system_trans_vectors(this->displacement_dh);
#else
  parallel::distributed::
    SolutionTransfer<dim, TrilinosWrappers::MPI::Vector, spacedim>
      system_trans_vectors(this->displacement_dh);
#endif

  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;
  displacement_relevant = displacement;
  sol_set_transfer.push_back(&displacement_relevant);

  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
        this->solid_tria.get()))
    {
      std::string triangulationName = prefix + ".triangulation";
      tria->save(prefix + ".triangulation");
    }
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::read_checkpoint(std::string prefix)
{
  // Setup an un-refined triangulation before loading
  setup_triangulation(true);

  // Read the triangulation from the checkpoint
  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string("You are trying to read a solid triangulation, "
                              "but the restart file <") +
                  filename + "> does not appear to exist!"));

  try
    {
      if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
            this->solid_tria.get()))
        tria->load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }

  // Setup dof-handler for solid and displacement
  solid_dh.distribute_dofs(*fe);
  setup_displacement();


  // Read displacement vector
  std::vector<TrilinosWrappers::MPI::Vector *> x_system(1);
  x_system[0] = &(displacement);

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::Vector,
                                          DoFHandler<dim, spacedim>>
    system_trans_vectors(this->displacement_dh);
#else
  parallel::distributed::
    SolutionTransfer<dim, TrilinosWrappers::MPI::Vector, spacedim>
      system_trans_vectors(this->displacement_dh);
#endif

  system_trans_vectors.deserialize(x_system);
  displacement_relevant = displacement;

  // Reset triangulation position using displacement vector
  move_solid_triangulation_with_displacement();


  // We did not checkpoint particles, we re-create them from scratch
  setup_particles();
}


template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::setup_displacement()
{
  displacement_dh.distribute_dofs(*this->displacement_fe);
  locally_owned_dofs = displacement_dh.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(displacement_dh,
                                          locally_relevant_dofs);

  displacement.reinit(locally_owned_dofs, mpi_communicator);
  displacement_relevant.reinit(locally_owned_dofs,
                               locally_relevant_dofs,
                               mpi_communicator);

  displacement = 0;
}


// Pre-compile the 2D, 3D and the 2D in 3D versions with the types that can
// occur
template class SolidBase<2>;
template class SolidBase<3>;
template class SolidBase<2, 3>;
