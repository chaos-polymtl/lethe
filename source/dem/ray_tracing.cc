// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/ray_tracing.h>

template <int dim, typename PropertiesIndex>
RayTracingSolver<dim, PropertiesIndex>::RayTracingSolver(
  RayTracingSolverParameters<dim> parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , parameters(parameters)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
  , photon_handler(triangulation, mapping, 4)
{}

template <int dim, typename PropertiesIndex>
void
RayTracingSolver<dim, PropertiesIndex>::setup_parameters()
{
  // Print simulation starting information
  pcout << std::endl;
  std::stringstream ss;
  ss << "Running on " << n_mpi_processes << " rank(s)";
  announce_string(pcout, ss.str(), '*');

  // Check if the output directory exists
  std::string output_dir_name = parameters.simulation_control.output_folder;
  struct stat buffer;

  // If output directory does not exist, create it
  if (this_mpi_process == 0)
    {
      if (stat(output_dir_name.c_str(), &buffer) != 0)
        {
          create_output_folder(output_dir_name);
        }
    }

  AdaptiveSparseContacts<dim, PropertiesIndex> dummy_sparse_contacts_object;

  // Set the simulation control as transient DEM
  simulation_control = std::make_shared<SimulationControlTransientDEM>(
    parameters.simulation_control);

  // Setup load balancing parameters and attach the correct functions to the
  // signals inside the triangulation
  load_balancing.set_parameters(parameters.model_parameters);
  load_balancing.copy_references(simulation_control,
                                 triangulation,
                                 particle_handler,
                                 dummy_sparse_contacts_object);
  load_balancing.connect_weight_signals();

  // Set the adaptive sparse contacts parameters

  // Set up the solid objects
  // setup_solid_objects();
}

template <int dim, typename PropertiesIndex>
std::shared_ptr<Insertion<dim, PropertiesIndex>>
RayTracingSolver<dim, PropertiesIndex>::set_particle_insertion_type()
{
  using namespace Parameters::Lagrangian;
  typename InsertionInfo<dim>::InsertionMethod insertion_method =
    parameters.insertion_info.insertion_method;

  std::vector<std::shared_ptr<Distribution>>
    dummy_size_distribution_object_container;

  switch (insertion_method)
    {
      case InsertionInfo<dim>::InsertionMethod::file:
        {
          return std::make_shared<InsertionFile<dim, PropertiesIndex>>(
            dummy_size_distribution_object_container,
            triangulation,
            parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::list:
        {
          return std::make_shared<InsertionList<dim, PropertiesIndex>>(
            dummy_size_distribution_object_container,
            triangulation,
            parameters);
        }
      default:
        throw(std::runtime_error("Invalid insertion method."));
    }
}

template <int dim, typename PropertiesIndex>
void
RayTracingSolver<dim, PropertiesIndex>::load_balance()
{
  TimerOutput::Scope t(this->computing_timer, "Load balancing");
  // Prepare particle handler for the adaptation of the triangulation to the
  // load
  particle_handler.prepare_for_coarsening_and_refinement();

  pcout << "-->Repartitionning triangulation" << std::endl;
  triangulation.repartition();

  // Unpack the particle handler after the mesh has been repartitioned
  particle_handler.unpack_after_coarsening_and_refinement();


  const auto average_minimum_maximum_cells =
    Utilities::MPI::min_max_avg(triangulation.n_active_cells(),
                                mpi_communicator);

  const auto average_minimum_maximum_particles =
    Utilities::MPI::min_max_avg(particle_handler.n_locally_owned_particles(),
                                mpi_communicator);

  pcout << "Load balance finished " << std::endl;
  pcout
    << "Average, minimum and maximum number of particles on the processors are "
    << average_minimum_maximum_particles.avg << " , "
    << average_minimum_maximum_particles.min << " and "
    << average_minimum_maximum_particles.max << std::endl;
  pcout << "Minimum and maximum number of cells owned by the processors are "
        << average_minimum_maximum_cells.min << " and "
        << average_minimum_maximum_cells.max << std::endl;

  setup_background_dofs();
}

template <int dim, typename PropertiesIndex>
void
RayTracingSolver<dim, PropertiesIndex>::insert_particles_and_photons()
{
  // Insert particles using the insertion object.
  particle_insertion_object->insert(particle_handler,
                                    triangulation,
                                    parameters);
}

template <int dim, typename PropertiesIndex>
void
RayTracingSolver<dim, PropertiesIndex>::finish_simulation()
{
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    this->computing_timer.print_summary();

  // Testing
  if (parameters.test.enabled)
    {
      // This needs to be coded later.
    }
}

template <int dim, typename PropertiesIndex>
void
RayTracingSolver<dim, PropertiesIndex>::solve()
{
  // Set up the parameters
  setup_parameters();
}
