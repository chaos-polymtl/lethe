/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Shahab Golshan, Polytechnique Montreal, 2019-
 */

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <core/solutions_output.h>
#include <dem/dem.h>

template <int dim>
DEMSolver<dim>::DEMSolver(DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
  , parameters(dem_parameters)
  , triangulation(this->mpi_communicator)
  , property_pool(DEM::get_number_properties())
  , mapping(1)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , particle_handler(triangulation, mapping, DEM::get_number_properties())
  , number_of_steps(parameters.simulation_control.timeEnd < 0.25 ? 20000 :
                                                                   50000)
{
  // Change the behavior of the timer for situations when you don't want outputs
  if (parameters.timer.type == Parameters::Timer::Type::none)
    computing_timer.disable_output();
}

// REFACTOR
// This function should be replaced by the function in the core/grids.h
template <int dim>
void
DEMSolver<dim>::read_mesh()
{
  // GMSH input
  if (parameters.mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(parameters.mesh.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (parameters.mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        triangulation,
        parameters.mesh.grid_type,
        parameters.mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  const int initial_size = parameters.mesh.initial_refinement;
  triangulation.refine_global(initial_size);
}

template <int dim>
void
DEMSolver<dim>::reinitialize_force(
  Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting properties of particle as local variable
      auto particle_properties = particle->get_properties();

      // Reinitializing forces and momentums of particles in the system
      particle_properties[DEM::PropertiesIndex::force_x] = 0;
      particle_properties[DEM::PropertiesIndex::force_y] = 0;

      particle_properties[DEM::PropertiesIndex::M_x] = 0;
      particle_properties[DEM::PropertiesIndex::M_y] = 0;

      if (dim == 3)
        {
          particle_properties[DEM::PropertiesIndex::force_z] = 0;
          particle_properties[DEM::PropertiesIndex::M_z]     = 0;
        }
    }
}

template <int dim>
std::shared_ptr<Insertion<dim>>
DEMSolver<dim>::set_insertion_type(const DEMSolverParameters<dim> &parameters)
{
  if (parameters.insertionInfo.insertion_method ==
      Parameters::Lagrangian::InsertionInfo::InsertionMethod::uniform)
    {
      insertion_object = std::make_shared<UniformInsertion<dim>>(parameters);
    }
  else if (parameters.insertionInfo.insertion_method ==
           Parameters::Lagrangian::InsertionInfo::InsertionMethod::non_uniform)
    {
      insertion_object = std::make_shared<NonUniformInsertion<dim>>(parameters);
    }
  else
    {
      throw "The chosen insertion method is invalid";
    }
  return insertion_object;
}

template <int dim>
std::shared_ptr<Integrator<dim>>
DEMSolver<dim>::set_integrator_type(const DEMSolverParameters<dim> &parameters)
{
  if (parameters.model_parmeters.integration_method ==
      Parameters::Lagrangian::ModelParameters::IntegrationMethod::
        velocity_verlet)
    {
      integrator_object = std::make_shared<VelocityVerletIntegrator<dim>>();
    }
  else if (parameters.model_parmeters.integration_method ==
           Parameters::Lagrangian::ModelParameters::IntegrationMethod::
             explicit_euler)
    {
      integrator_object = std::make_shared<ExplicitEulerIntegrator<dim>>();
    }
  else
    {
      throw "The chosen integration method is invalid";
    }
  return integrator_object;
}

template <int dim>
std::shared_ptr<PPContactForce<dim>>
DEMSolver<dim>::set_pp_contact_force(const DEMSolverParameters<dim> &parameters)
{
  if (parameters.model_parmeters.pp_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PPContactForceModel::pp_linear)
    {
      pp_contact_force_object = std::make_shared<PPLinearForce<dim>>();
    }
  else if (parameters.model_parmeters.pp_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::PPContactForceModel::
             pp_nonlinear)
    {
      pp_contact_force_object = std::make_shared<PPNonLinearForce<dim>>();
    }
  else
    {
      throw "The chosen particle-particle contact force model is invalid";
    }
  return pp_contact_force_object;
}

template <int dim>
std::shared_ptr<PWContactForce<dim>>
DEMSolver<dim>::set_pw_contact_force(const DEMSolverParameters<dim> &parameters)
{
  if (parameters.model_parmeters.pw_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PWContactForceModel::pw_linear)
    {
      pw_contact_force_object = std::make_shared<PWLinearForce<dim>>();
    }
  else if (parameters.model_parmeters.pw_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::PWContactForceModel::
             pw_nonlinear)
    {
      pw_contact_force_object = std::make_shared<PWNonLinearForce<dim>>();
    }
  else
    {
      throw "The chosen particle-wall contact force model is invalid";
    }
  return pw_contact_force_object;
}

template <int dim>
std::map<int, Particles::ParticleIterator<dim>>
DEMSolver<dim>::update_particle_container(
  const Particles::ParticleHandler<dim> *particle_handler)
{
  std::map<int, Particles::ParticleIterator<dim>> particle_container;
  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  return particle_container;
}

template <int dim>
void
DEMSolver<dim>::update_pp_contact_container_iterators(
  std::map<int, std::map<int, pp_contact_info_struct<dim>>> &adjacent_particles,
  const std::map<int, Particles::ParticleIterator<dim>> &    particle_container)
{
  for (auto adjacent_particles_iterator = adjacent_particles.begin();
       adjacent_particles_iterator != adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();
           ++pp_map_iterator)
        {
          int particle_two_id = pp_map_iterator->first;
          pp_map_iterator->second.particle_one =
            particle_container.at(particle_one_id);
          pp_map_iterator->second.particle_two =
            particle_container.at(particle_two_id);
        }
    }
}

template <int dim>
void
DEMSolver<dim>::update_pw_contact_container_iterators(
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    &                                                    pw_pairs_in_contact,
  const std::map<int, Particles::ParticleIterator<dim>> &particle_container)
{
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator)
    {
      int particle_id = pw_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content = &pw_pairs_in_contact_iterator->second;

      for (auto pw_map_iterator = pairs_in_contant_content->begin();
           pw_map_iterator != pairs_in_contant_content->end();
           ++pw_map_iterator)
        {
          pw_map_iterator->second.particle = particle_container.at(particle_id);
        }
    }
}

template <int dim>
void
DEMSolver<dim>::update_particle_point_line_contact_container_iterators(
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_points_in_contact,
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_lines_in_contact,
  const std::map<int, Particles::ParticleIterator<dim>> &particle_container)
{
  for (auto particle_point_pairs_in_contact_iterator =
         particle_points_in_contact.begin();
       particle_point_pairs_in_contact_iterator !=
       particle_points_in_contact.end();
       ++particle_point_pairs_in_contact_iterator)
    {
      int  particle_id = particle_point_pairs_in_contact_iterator->first;
      auto pairs_in_contant_content =
        &particle_point_pairs_in_contact_iterator->second;
      pairs_in_contant_content->particle = particle_container.at(particle_id);
    }

  for (auto particle_line_pairs_in_contact_iterator =
         particle_lines_in_contact.begin();
       particle_line_pairs_in_contact_iterator !=
       particle_lines_in_contact.end();
       ++particle_line_pairs_in_contact_iterator)
    {
      int  particle_id = particle_line_pairs_in_contact_iterator->first;
      auto pairs_in_contant_content =
        &particle_line_pairs_in_contact_iterator->second;
      pairs_in_contant_content->particle = particle_container.at(particle_id);
    }
}

template <int dim>
void
DEMSolver<dim>::write_output_results()
{
  const std::string  folder = parameters.simulation_control.output_folder;
  const std::string  solution_name = parameters.simulation_control.output_name;
  const unsigned int iter          = DEM_step;
  const double       time          = DEM_time;
  const unsigned int group_files   = parameters.simulation_control.group_files;

  Visualization<dim> particle_data_out;
  particle_data_out.build_patches(particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(pvdhandler,
                            particle_data_out,
                            folder,
                            solution_name,
                            time,
                            iter,
                            group_files,
                            mpi_communicator);
}

template <int dim>
void
DEMSolver<dim>::solve()
{
  // Reading mesh
  read_mesh();

  // Initializing variables
  Tensor<1, dim> g;

  g[0] = parameters.physicalProperties.gx;
  g[1] = parameters.physicalProperties.gy;
  if (dim == 3)
    {
      g[2] = parameters.physicalProperties.gz;
    }

  auto properties = properties_class.get_properties_name();

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbors_object;
  cell_neighbor_list = cell_neighbors_object.find_cell_neighbors(triangulation);

  // Finding boundary cells with faces
  FindBoundaryCellsInformation<dim> boundary_cell_object;
  boundary_cells_information =
    boundary_cell_object.find_boundary_cells_information(
      boundary_cells_with_faces, triangulation);

  // Finding boundary cells with lines and points
  boundary_cell_object.find_particle_point_and_line_contact_cells(
    boundary_cells_with_faces,
    triangulation,
    boundary_cells_with_lines,
    boundary_cells_with_points);

  // Setting chosen contact force, insertion and integration methods
  insertion_object        = set_insertion_type(parameters);
  integrator_object       = set_integrator_type(parameters);
  pp_contact_force_object = set_pp_contact_force(parameters);
  pw_contact_force_object = set_pw_contact_force(parameters);

  // DEM engine iterator:
  while (DEM_step < number_of_steps)
    {
      // Moving walls

      // Insertion
      // Defining a bool variable to specify the insertion steps
      bool insertion_step = 0;

      // REFACTORING
      // Why is this a floating point modulu? Should be a regular modulo for
      // integers... It would be better if the insertion class would control
      // these elements You could have a mother class for insertion and derive
      // from there the type of insertion
      if (fmod(DEM_step, parameters.insertionInfo.insertion_frequency) == 1)
        {
          computing_timer.enter_subsection("insertion");
          // REFACTORING
          insertion_object->insert(particle_handler, triangulation, parameters);
          insertion_step = 1;
          computing_timer.leave_subsection();
        }

      // Sort particles in cells
      if (insertion_step == 1 || DEM_step % pp_broad_search_frequency == 0 ||
          DEM_step % pw_broad_search_frequency == 0)
        {
          computing_timer.enter_subsection("sort_particles_in_cells");
          particle_handler.sort_particles_into_subdomains_and_cells();
          particle_container.clear();
          particle_container = update_particle_container(&particle_handler);
          update_pp_contact_container_iterators(adjacent_particles,
                                                particle_container);
          update_pw_contact_container_iterators(pw_pairs_in_contact,
                                                particle_container);
          update_particle_point_line_contact_container_iterators(
            particle_points_in_contact,
            particle_lines_in_contact,
            particle_container);
          computing_timer.leave_subsection();
        }

      // Force reinitilization
      computing_timer.enter_subsection("reinitialize_forces");
      reinitialize_force(particle_handler);
      computing_timer.leave_subsection();

      // PP contact search
      // PP broad search
      if (insertion_step == 1 || DEM_step % pp_broad_search_frequency == 0)
        {
          computing_timer.enter_subsection("pp_broad_search");
          pp_broad_search_object.find_PP_Contact_Pairs(particle_handler,
                                                       cell_neighbor_list,
                                                       contact_pair_candidates);
          computing_timer.leave_subsection();
        }

      // PP fine search
      if (insertion_step == 1 || DEM_step % pp_fine_search_frequency == 0)
        {
          computing_timer.enter_subsection("pp_fine_search");
          pp_fine_search_object.pp_Fine_Search(contact_pair_candidates,
                                               adjacent_particles,
                                               neighborhood_threshold);
          computing_timer.leave_subsection();
        }

      // PP contact force
      computing_timer.enter_subsection("pp_contact_force");
      pp_contact_force_object->calculate_pp_contact_force(&adjacent_particles,
                                                          parameters,
                                                          dt);
      computing_timer.leave_subsection();

      // PW contact search
      // PW broad search
      if (insertion_step == 1 || DEM_step % pw_broad_search_frequency == 0)
        {
          computing_timer.enter_subsection("pw_broad_search");
          pw_broad_search_object.find_PW_Contact_Pairs(
            boundary_cells_information,
            particle_handler,
            pw_contact_candidates);

          particle_point_contact_candidates =
            particle_point_line_broad_search_object
              .find_Particle_Point_Contact_Pairs(particle_handler,
                                                 boundary_cells_with_points);
          if (dim == 3)
            {
              particle_line_contact_candidates =
                particle_point_line_broad_search_object
                  .find_Particle_Line_Contact_Pairs(particle_handler,
                                                    boundary_cells_with_lines);
            }
          computing_timer.leave_subsection();
        }

      // PW fine search
      computing_timer.enter_subsection("pw_fine_search");
      pw_fine_search_object.pw_Fine_Search(pw_contact_candidates,
                                           pw_pairs_in_contact,
                                           dt);
      particle_points_in_contact =
        particle_point_line_fine_search_object.Particle_Point_Fine_Search(
          particle_point_contact_candidates);
      if (dim == 3)
        {
          particle_lines_in_contact =
            particle_point_line_fine_search_object.Particle_Line_Fine_Search(
              particle_line_contact_candidates);
        }
      computing_timer.leave_subsection();

      // PW contact force:
      computing_timer.enter_subsection("pw_contact_force");
      pw_contact_force_object->calculate_pw_contact_force(&pw_pairs_in_contact,
                                                          parameters);
      particle_point_line_contact_force_object
        .calculate_particle_point_line_contact_force(
          &particle_points_in_contact, parameters);
      if (dim == 3)
        {
          particle_point_line_contact_force_object
            .calculate_particle_point_line_contact_force(
              &particle_lines_in_contact, parameters);
        }
      computing_timer.leave_subsection();

      // Integration
      computing_timer.enter_subsection("integration");
      integrator_object->integrate(particle_handler, g, dt);
      computing_timer.leave_subsection();

      // Visualization
      computing_timer.enter_subsection("visualization");
      // REFACTORING
      // Should be put inside a function
      if (DEM_step % parameters.simulation_control.output_frequency == 0)
        {
          write_output_results();
        }
      computing_timer.leave_subsection();

      // Print iteration
      // REFACTORING
      // Should be put into a simulation control type of object
      if (DEM_step % print_info_frequency == 0)
        {
          std::cout << "Step " << DEM_step << std::endl;
        }

      // Update:
      DEM_step = DEM_step + 1;
      DEM_time = DEM_step * dt;
    }

  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    {
      this->computing_timer.print_summary();
    }

  // Testing
  if (parameters.test.enabled)
    {
      visualization_object.print_xyz(particle_handler, properties);
    }

  /*
    while (parameters.simulation_control.integrate()) {
      printTime(this->pcout, parameters.simulation_control);

    }
        */
}

template class DEMSolver<2>;
template class DEMSolver<3>;
