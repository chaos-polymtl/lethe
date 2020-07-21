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
    : mpi_communicator(MPI_COMM_WORLD),
      n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)),
      this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)),
      pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0}),
      parameters(dem_parameters), triangulation(this->mpi_communicator),
      property_pool(DEM::get_number_properties()), mapping(1),
      computing_timer(this->mpi_communicator, this->pcout, TimerOutput::summary,
                      TimerOutput::wall_times),
      particle_handler(triangulation, mapping, DEM::get_number_properties()),
      background_dh(triangulation) {
  // Change the behavior of the timer for situations when you don't want outputs
  if (parameters.timer.type == Parameters::Timer::Type::none)
    computing_timer.disable_output();

  simulation_control = std::make_shared<SimulationControlTransientDEM>(
      parameters.simulation_control);
}

template <int dim> void DEMSolver<dim>::read_mesh() {
  // GMSH input
  if (parameters.mesh.type == Parameters::Mesh::Type::gmsh) {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file(parameters.mesh.file_name);
    grid_in.read_msh(input_file);
  }

  // Dealii grids
  else if (parameters.mesh.type == Parameters::Mesh::Type::dealii) {
    GridGenerator::generate_from_name_and_arguments(
        triangulation, parameters.mesh.grid_type,
        parameters.mesh.grid_arguments);
  } else
    throw std::runtime_error(
        "Unsupported mesh type - mesh will not be created");

  const int initial_size = parameters.mesh.initial_refinement;
  triangulation.refine_global(initial_size);
}

template <int dim> void DEMSolver<dim>::setup_background_dofs() {
  FE_Q<dim> background_fe(1);
  background_dh.distribute_dofs(background_fe);
}

template <int dim> bool DEMSolver<dim>::insert_particles() {
  TimerOutput::Scope t(computing_timer, "insertion");

  if (fmod(simulation_control->get_step_number(),
           parameters.insertion_info.insertion_frequency) == 1) {
    insertion_object->insert(particle_handler, triangulation, parameters);
    return true;
  }
  return false;
}

template <int dim> void DEMSolver<dim>::clear_contact_containers() {

  // std::cout << "size " << local_adjacent_particles.size() << std::endl;

  cleared_local_adjacent_particles.clear();
  cleared_ghost_adjacent_particles.clear();
  cleared_pw_pairs_in_contact.clear();

  for (auto adjacent_particles_iterator = local_adjacent_particles.begin();
       adjacent_particles_iterator != local_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    int particle_one_id = adjacent_particles_iterator->first;

    auto pairs_in_contant_content = &adjacent_particles_iterator->second;
    for (auto pp_map_iterator = pairs_in_contant_content->begin();
         pp_map_iterator != pairs_in_contant_content->end();
         ++pp_map_iterator) {

      int particle_two_id = pp_map_iterator->first;
      // std::cout << "size " << local_contact_pair_candidates.size() <<
      // std::endl;

      //  std::cout << pp_map_iterator->second.tangential_overlap << std::endl;

      auto search_iterator_one = local_contact_pair_candidates.find(
          std::make_pair(particle_one_id, particle_two_id));
      auto search_iterator_two = local_contact_pair_candidates.find(
          std::make_pair(particle_two_id, particle_one_id));

      if (search_iterator_one != local_contact_pair_candidates.end()) {
        // std::cout << pp_map_iterator->second.tangential_overlap << std::endl;
        local_contact_pair_candidates.erase(search_iterator_one);
        cleared_local_adjacent_particles[particle_one_id].insert(
            {particle_two_id, pp_map_iterator->second});

      } else if (search_iterator_two != local_contact_pair_candidates.end()) {
        //   std::cout << "in second if " << std::endl;
        local_contact_pair_candidates.erase(search_iterator_two);
        cleared_local_adjacent_particles[particle_one_id].insert(
            {particle_two_id, pp_map_iterator->second});
      }
    }
  }

  // The same for local-ghost particle containers
  for (auto adjacent_particles_iterator = ghost_adjacent_particles.begin();
       adjacent_particles_iterator != ghost_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    int particle_one_id = adjacent_particles_iterator->first;

    auto pairs_in_contant_content = &adjacent_particles_iterator->second;
    for (auto pp_map_iterator = pairs_in_contant_content->begin();
         pp_map_iterator != pairs_in_contant_content->end();
         ++pp_map_iterator) {

      int particle_two_id = pp_map_iterator->first;

      auto search_iterator_one = ghost_contact_pair_candidates.find(
          std::make_pair(particle_one_id, particle_two_id));
      auto search_iterator_two = ghost_contact_pair_candidates.find(
          std::make_pair(particle_two_id, particle_one_id));

      if (search_iterator_one != ghost_contact_pair_candidates.end()) {
        ghost_contact_pair_candidates.erase(search_iterator_one);
        cleared_ghost_adjacent_particles[particle_one_id].insert(
            {particle_two_id, pp_map_iterator->second});

      } else if (search_iterator_two != ghost_contact_pair_candidates.end()) {
        ghost_contact_pair_candidates.erase(search_iterator_two);
        cleared_ghost_adjacent_particles[particle_one_id].insert(
            {particle_two_id, pp_map_iterator->second});
      }
    }
  }

  // Particle-wall contacts
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator) {

    int particle_id = pw_pairs_in_contact_iterator->first;

    auto pairs_in_contant_content = &pw_pairs_in_contact_iterator->second;

    for (auto pw_map_iterator = pairs_in_contant_content->begin();
         pw_map_iterator != pairs_in_contant_content->end();
         ++pw_map_iterator) {

      int face_id = pw_map_iterator->first;

      auto search_iterator =
          pw_contact_candidates.find(std::make_pair(particle_id, face_id));

      if (search_iterator != pw_contact_candidates.end()) {
        pw_contact_candidates.erase(search_iterator);
        cleared_pw_pairs_in_contact[particle_id].insert(
            {face_id, pw_map_iterator->second});
      }
    }
  }
}

template <int dim> void DEMSolver<dim>::locate_particles_in_cells() {
  computing_timer.enter_subsection("sort_particles_in_cells");

  particle_container.clear();
  update_particle_container(particle_container, &particle_handler);
  // std::cout << "after container " << std::endl;
  update_pp_contact_container_iterators(cleared_local_adjacent_particles,
                                        cleared_ghost_adjacent_particles,
                                        particle_container);
  // std::cout << "after update pp " << std::endl;
  update_pw_contact_container_iterators(cleared_pw_pairs_in_contact,
                                        particle_container);
  //  std::cout << "after update pw " << std::endl;

  update_particle_point_line_contact_container_iterators(
      particle_points_in_contact, particle_lines_in_contact,
      particle_container);

  computing_timer.leave_subsection();
}

template <int dim> void DEMSolver<dim>::particle_wall_broad_search() {
  computing_timer.enter_subsection("pw_broad_search");
  pw_broad_search_object.find_PW_Contact_Pairs(
      boundary_cells_information, particle_handler, pw_contact_candidates);

  particle_point_contact_candidates =
      particle_point_line_broad_search_object.find_Particle_Point_Contact_Pairs(
          particle_handler, boundary_cells_with_points);
  if (dim == 3) {
    particle_line_contact_candidates =
        particle_point_line_broad_search_object
            .find_Particle_Line_Contact_Pairs(particle_handler,
                                              boundary_cells_with_lines);
  }
  computing_timer.leave_subsection();
}

template <int dim> void DEMSolver<dim>::particle_wall_fine_search() {
  computing_timer.enter_subsection("pw_fine_search");
  pw_fine_search_object.pw_Fine_Search(pw_contact_candidates,
                                       cleared_pw_pairs_in_contact,
                                       simulation_control->get_time_step());
  particle_points_in_contact =
      particle_point_line_fine_search_object.Particle_Point_Fine_Search(
          particle_point_contact_candidates);
  if (dim == 3) {
    particle_lines_in_contact =
        particle_point_line_fine_search_object.Particle_Line_Fine_Search(
            particle_line_contact_candidates);
  }
  computing_timer.leave_subsection();
}

template <int dim> void DEMSolver<dim>::particle_wall_contact_force() {
  computing_timer.enter_subsection("pw_contact_force");
  pw_contact_force_object->calculate_pw_contact_force(
      &cleared_pw_pairs_in_contact, parameters);
  particle_point_line_contact_force_object
      .calculate_particle_point_line_contact_force(&particle_points_in_contact,
                                                   parameters);

  if (dim == 3) {
    particle_point_line_contact_force_object
        .calculate_particle_point_line_contact_force(&particle_lines_in_contact,
                                                     parameters);
  }
  computing_timer.leave_subsection();
}

template <int dim> void DEMSolver<dim>::finish_simulation() {
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    this->computing_timer.print_summary();

  // Testing
  if (parameters.test.enabled) {
    auto properties = properties_class.get_properties_name();
    Visualization<dim> visualization_object;
    visualization_object.print_xyz(particle_handler, properties);
  }
}

template <int dim>
void DEMSolver<dim>::reinitialize_force(
    Particles::ParticleHandler<dim> &particle_handler) {
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    // Getting properties of particle as local variable
    auto particle_properties = particle->get_properties();

    // Reinitializing forces and momentums of particles in the system
    particle_properties[DEM::PropertiesIndex::force_x] = 0;
    particle_properties[DEM::PropertiesIndex::force_y] = 0;

    particle_properties[DEM::PropertiesIndex::M_x] = 0;
    particle_properties[DEM::PropertiesIndex::M_y] = 0;

    if (dim == 3) {
      particle_properties[DEM::PropertiesIndex::force_z] = 0;
      particle_properties[DEM::PropertiesIndex::M_z] = 0;
    }
  }
}

template <int dim>
std::shared_ptr<Insertion<dim>>
DEMSolver<dim>::set_insertion_type(const DEMSolverParameters<dim> &parameters) {
  if (parameters.insertion_info.insertion_method ==
      Parameters::Lagrangian::InsertionInfo::InsertionMethod::uniform) {
    insertion_object = std::make_shared<UniformInsertion<dim>>(parameters);
  } else if (parameters.insertion_info.insertion_method ==
             Parameters::Lagrangian::InsertionInfo::InsertionMethod::
                 non_uniform) {
    insertion_object = std::make_shared<NonUniformInsertion<dim>>(parameters);
  } else {
    throw "The chosen insertion method is invalid";
  }
  return insertion_object;
}

template <int dim>
std::shared_ptr<Integrator<dim>> DEMSolver<dim>::set_integrator_type(
    const DEMSolverParameters<dim> &parameters) {
  if (parameters.model_parameters.integration_method ==
      Parameters::Lagrangian::ModelParameters::IntegrationMethod::
          velocity_verlet) {
    integrator_object = std::make_shared<VelocityVerletIntegrator<dim>>();
  } else if (parameters.model_parameters.integration_method ==
             Parameters::Lagrangian::ModelParameters::IntegrationMethod::
                 explicit_euler) {
    integrator_object = std::make_shared<ExplicitEulerIntegrator<dim>>();
  } else {
    throw "The chosen integration method is invalid";
  }
  return integrator_object;
}

template <int dim>
std::shared_ptr<PPContactForce<dim>> DEMSolver<dim>::set_pp_contact_force(
    const DEMSolverParameters<dim> &parameters) {
  if (parameters.model_parameters.pp_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PPContactForceModel::pp_linear) {
    pp_contact_force_object = std::make_shared<PPLinearForce<dim>>();
  } else if (parameters.model_parameters.pp_contact_force_method ==
             Parameters::Lagrangian::ModelParameters::PPContactForceModel::
                 pp_nonlinear) {
    pp_contact_force_object = std::make_shared<PPNonLinearForce<dim>>();
  } else {
    throw "The chosen particle-particle contact force model is invalid";
  }
  return pp_contact_force_object;
}

template <int dim>
std::shared_ptr<PWContactForce<dim>> DEMSolver<dim>::set_pw_contact_force(
    const DEMSolverParameters<dim> &parameters) {
  if (parameters.model_parameters.pw_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PWContactForceModel::pw_linear) {
    pw_contact_force_object = std::make_shared<PWLinearForce<dim>>();
  } else if (parameters.model_parameters.pw_contact_force_method ==
             Parameters::Lagrangian::ModelParameters::PWContactForceModel::
                 pw_nonlinear) {
    pw_contact_force_object = std::make_shared<PWNonLinearForce<dim>>();
  } else {
    throw "The chosen particle-wall contact force model is invalid";
  }
  return pw_contact_force_object;
}

template <int dim>
void DEMSolver<dim>::update_particle_container(
    std::map<int, Particles::ParticleIterator<dim>> &particle_container,
    Particles::ParticleHandler<dim> *particle_handler) {

  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end(); ++particle_iterator) {
    particle_container[particle_iterator->get_id()] = particle_iterator;
  }

  // Loop over the ghost particles
  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator) {
    particle_container[particle_iterator->get_id()] = particle_iterator;
  }
}

template <int dim>
void DEMSolver<dim>::update_pp_contact_container_iterators(
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
        &cleared_local_adjacent_particles,
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
        &cleared_ghost_adjacent_particles,
    const std::map<int, Particles::ParticleIterator<dim>> &particle_container) {

  for (auto adjacent_particles_iterator =
           cleared_local_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_local_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    int particle_one_id = adjacent_particles_iterator->first;
    auto pairs_in_contant_content = &adjacent_particles_iterator->second;
    for (auto pp_map_iterator = pairs_in_contant_content->begin();
         pp_map_iterator != pairs_in_contant_content->end();
         ++pp_map_iterator) {
      int particle_two_id = pp_map_iterator->first;

      pp_map_iterator->second.particle_one =
          particle_container.at(particle_one_id);
      pp_map_iterator->second.particle_two =
          particle_container.at(particle_two_id);
    }
  }

  // Doing the same thing for ghost container
  for (auto adjacent_particles_iterator =
           cleared_ghost_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_ghost_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    int particle_one_id = adjacent_particles_iterator->first;
    auto pairs_in_contant_content = &adjacent_particles_iterator->second;
    for (auto pp_map_iterator = pairs_in_contant_content->begin();
         pp_map_iterator != pairs_in_contant_content->end();
         ++pp_map_iterator) {
      int particle_two_id = pp_map_iterator->first;

      pp_map_iterator->second.particle_one =
          particle_container.at(particle_one_id);
      pp_map_iterator->second.particle_two =
          particle_container.at(particle_two_id);
    }
  }
}

template <int dim>
void DEMSolver<dim>::update_pw_contact_container_iterators(
    std::map<int, std::map<int, pw_contact_info_struct<dim>>>
        &cleared_pw_pairs_in_contact,
    const std::map<int, Particles::ParticleIterator<dim>> &particle_container) {

  for (auto pw_pairs_in_contact_iterator = cleared_pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != cleared_pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator) {

    int particle_id = pw_pairs_in_contact_iterator->first;

    auto pairs_in_contant_content = &pw_pairs_in_contact_iterator->second;

    for (auto pw_map_iterator = pairs_in_contant_content->begin();
         pw_map_iterator != pairs_in_contant_content->end();
         ++pw_map_iterator) {

      pw_map_iterator->second.particle = particle_container.at(particle_id);
    }
  }
}

template <int dim>
void DEMSolver<dim>::update_particle_point_line_contact_container_iterators(
    std::map<int, particle_point_line_contact_info_struct<dim>>
        &particle_points_in_contact,
    std::map<int, particle_point_line_contact_info_struct<dim>>
        &particle_lines_in_contact,
    const std::map<int, Particles::ParticleIterator<dim>> &particle_container) {
  for (auto particle_point_pairs_in_contact_iterator =
           particle_points_in_contact.begin();
       particle_point_pairs_in_contact_iterator !=
       particle_points_in_contact.end();
       ++particle_point_pairs_in_contact_iterator) {
    int particle_id = particle_point_pairs_in_contact_iterator->first;
    auto pairs_in_contact_content =
        &particle_point_pairs_in_contact_iterator->second;
    pairs_in_contact_content->particle = particle_container.at(particle_id);
  }

  for (auto particle_line_pairs_in_contact_iterator =
           particle_lines_in_contact.begin();
       particle_line_pairs_in_contact_iterator !=
       particle_lines_in_contact.end();
       ++particle_line_pairs_in_contact_iterator) {
    int particle_id = particle_line_pairs_in_contact_iterator->first;
    auto pairs_in_contant_content =
        &particle_line_pairs_in_contact_iterator->second;
    pairs_in_contant_content->particle = particle_container.at(particle_id);
  }
}

template <int dim> void DEMSolver<dim>::write_output_results() {
  const std::string folder = parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
      parameters.simulation_control.output_name;
  const unsigned int iter = simulation_control->get_step_number();
  const double time = simulation_control->get_current_time();
  const unsigned int group_files = parameters.simulation_control.group_files;

  // Write particles
  Visualization<dim> particle_data_out;
  particle_data_out.build_patches(particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler, particle_data_out, folder,
                            particles_solution_name, time, iter, group_files,
                            mpi_communicator);

  // Write background grid
  DataOut<dim> background_data_out;

  background_data_out.attach_dof_handler(background_dh);

  // Attach the solution data to data_out object
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  background_data_out.add_data_vector(subdomain, "subdomain");

  const std::string grid_solution_name =
      parameters.simulation_control.output_name + "-grid";

  background_data_out.build_patches();

  write_vtu_and_pvd<dim>(grid_pvdhandler, background_data_out, folder,
                         grid_solution_name, time, iter, group_files,
                         mpi_communicator);
}

template <int dim> void DEMSolver<dim>::solve() {

  // Reading mesh
  read_mesh();
  //****************//setup_background_dofs();

  // Initialize DEM body force
  Tensor<1, dim> g;

  g[0] = parameters.physical_properties.gx;
  g[1] = parameters.physical_properties.gy;
  if (dim == 3) {
    g[2] = parameters.physical_properties.gz;
  }

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbors_object;
  cell_neighbors_object.find_cell_neighbors(
      triangulation, cells_local_neighbor_list, cells_ghost_neighbor_list);

  /*
 if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1) {
   for (auto iterator = cells_local_neighbor_list.begin();
        iterator != cells_local_neighbor_list.end(); ++iterator) {

     std::vector<typename Triangulation<dim>::active_cell_iterator> content =
         *iterator;
     typename Triangulation<dim>::active_cell_iterator first_cell =
         *content.begin();

     std::cout << "The main cell is " << first_cell->id() << " are: ";

     for (auto iterator = 0; iterator != content.size(); ++iterator) {
       std::cout << content[iterator]->id() << " , ";
     }
     std::cout << std::endl;
   }
 }

 if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1) {
   for (auto iterator = cells_ghost_neighbor_list.begin();
        iterator != cells_ghost_neighbor_list.end(); ++iterator) {

     std::vector<typename Triangulation<dim>::active_cell_iterator> content =
         *iterator;
     auto first_cell = *content.begin();
     std::cout << "(Ghost list) The main cell is " << first_cell->id()
               << " are: ";

     for (auto iterator = 0; iterator != content.size(); ++iterator) {
       std::cout << content[iterator]->id() << " , ";
     }
     std::cout << std::endl;
   }
 }
   */

  // Finding boundary cells with faces
  FindBoundaryCellsInformation<dim> boundary_cell_object;
  boundary_cells_information =
      boundary_cell_object.find_boundary_cells_information(
          boundary_cells_with_faces, triangulation);

  // Finding boundary cells with lines and points
  boundary_cell_object.find_particle_point_and_line_contact_cells(
      boundary_cells_with_faces, triangulation, boundary_cells_with_lines,
      boundary_cells_with_points);

  // Setting chosen contact force, insertion and integration methods
  insertion_object = set_insertion_type(parameters);
  integrator_object = set_integrator_type(parameters);
  pp_contact_force_object = set_pp_contact_force(parameters);
  pw_contact_force_object = set_pw_contact_force(parameters);

  const unsigned int pp_broad_search_frequency =
      parameters.model_parameters.pp_broad_search_frequency;
  const unsigned int pw_broad_search_frequency =
      parameters.model_parameters.pw_broad_search_frequency;
  const unsigned int pp_fine_search_frequency =
      parameters.model_parameters.pp_fine_search_frequency;

  // DEM engine iterator:
  while (simulation_control->integrate()) {

    simulation_control->print_progression(pcout);
    const unsigned int step_number = simulation_control->get_step_number();

    // Keep track if particles were inserted this step
    bool particles_were_inserted = insert_particles();

    // Sort particles in cells
    if (particles_were_inserted ||
        step_number % pp_broad_search_frequency == 0 ||
        step_number % pw_broad_search_frequency == 0) {
      particle_handler.sort_particles_into_subdomains_and_cells();
    }

    //  if (particles_were_inserted ||
    //      step_number % pp_broad_search_frequency == 0 ||
    //      step_number % pw_broad_search_frequency == 0) {
    particle_handler.exchange_ghost_particles();
    //   }

    //  pcout << "-> before force re-init" << std::endl;

    // Force reinitilization
    computing_timer.enter_subsection("reinitialize_forces");
    reinitialize_force(particle_handler);
    computing_timer.leave_subsection();

    //   pcout << "-> force re-init" << std::endl;

    //  pcout << "-> before pp broad" << std::endl;
    // Broad particle-particle contact search
    if (particles_were_inserted ||
        step_number % pp_broad_search_frequency == 0) {
      computing_timer.enter_subsection("pp_broad_search");
      pp_broad_search_object.find_PP_Contact_Pairs(
          particle_handler, &cells_local_neighbor_list,
          &cells_ghost_neighbor_list, local_contact_pair_candidates,
          ghost_contact_pair_candidates);
      computing_timer.leave_subsection();
    }

    // Particle-wall broad contact search
    if (particles_were_inserted || step_number % pw_broad_search_frequency == 0)
      particle_wall_broad_search();
    // pcout << "-> w-broad" << std::endl;

    //  if (particles_were_inserted ||
    //      step_number % pp_broad_search_frequency == 0 ||
    //      step_number % pw_broad_search_frequency == 0) {
    clear_contact_containers();
    //  }

    //  pcout << "-> before " << std::endl;

    //   if (particles_were_inserted ||
    //      step_number % pp_broad_search_frequency == 0 ||
    //       step_number % pw_broad_search_frequency == 0) {
    locate_particles_in_cells();
    //   }
    //  pcout << "-> cell located" << std::endl;

    //   pcout << "-> before pp fine" << std::endl;

    // Particle-particle fine search
    if (particles_were_inserted ||
        step_number % pp_fine_search_frequency == 0) {
      computing_timer.enter_subsection("pp_fine_search");
      const double neighborhood_threshold =
          parameters.model_parameters.neighborhood_threshold *
          parameters.physical_properties.diameter;
      pp_fine_search_object.pp_Fine_Search(
          local_contact_pair_candidates, ghost_contact_pair_candidates,
          cleared_local_adjacent_particles, cleared_ghost_adjacent_particles,
          neighborhood_threshold);
      computing_timer.leave_subsection();
    }

    // /*
    //  pcout << "-> before pp force" << std::endl;

    // Particle-particle contact force
    computing_timer.enter_subsection("pp_contact_force");
    pp_contact_force_object->calculate_pp_contact_force(
        &cleared_local_adjacent_particles, &cleared_ghost_adjacent_particles,
        parameters, simulation_control->get_time_step());
    computing_timer.leave_subsection();

    //  pcout << "-> p-force" << std::endl;

    // */

    // Particles-wall fine search
    particle_wall_fine_search();

    //   pcout << "-> w-fine" << std::endl;

    // Particles-walls contact force:
    particle_wall_contact_force();

    //  pcout << "-> w-force" << std::endl;

    // Integration
    computing_timer.enter_subsection("integration");
    integrator_object->integrate(particle_handler, g,
                                 simulation_control->get_time_step());
    computing_timer.leave_subsection();

    //  pcout << "-> integration" << std::endl;

    //      // Visualization
    if (simulation_control->is_output_iteration()) {
      computing_timer.enter_subsection("visualization");
      write_output_results();
      computing_timer.leave_subsection();
    }

    //   pcout << "-> viz" << std::endl;

    {
      local_adjacent_particles.clear();
      ghost_adjacent_particles.clear();
      cleared_pw_pairs_in_contact.clear();
      local_adjacent_particles = cleared_local_adjacent_particles;
      ghost_adjacent_particles = cleared_ghost_adjacent_particles;
      pw_pairs_in_contact = cleared_pw_pairs_in_contact;
    }
  }

  finish_simulation();
}

template class DEMSolver<2>;
template class DEMSolver<3>;
