#include <core/manifolds.h>
#include <core/solutions_output.h>

#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/distributions.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/force_chains_visualization.h>
#include <dem/gear3_integrator.h>
#include <dem/input_parameter_inspection.h>
#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>
#include <dem/post_processing.h>
#include <dem/read_checkpoint.h>
#include <dem/read_mesh.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/write_checkpoint.h>

#include <deal.II/base/table_handler.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_out.h>

#include <sys/stat.h>

#include <sstream>

template <int dim>
DEMSolver<dim>::DEMSolver(DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , parameters(dem_parameters)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , particles_insertion_step(0)
  , contact_build_number(0)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , particle_handler(triangulation, mapping, DEM::get_number_properties())
  , contact_detection_step(true)
  , load_balance_step(true)
  , checkpoint_step(true)
  , contact_detection_frequency(
      parameters.model_parameters.contact_detection_frequency)
  , insertion_frequency(parameters.insertion_info.insertion_frequency)
  , has_periodic_boundaries(false)
  , background_dh(triangulation)
  , has_solid_objects(false)
  , size_distribution_object_container(
      parameters.lagrangian_physical_properties.particle_type_number)
  , has_sparse_contacts(false)
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

  // Change the behavior of the timer for situations when you don't want outputs
  if (parameters.timer.type == Parameters::Timer::Type::none)
    computing_timer.disable_output();
  simulation_control = std::make_shared<SimulationControlTransientDEM>(
    parameters.simulation_control);

  // In order to consider the particles when repartitioning the triangulation
  // the algorithm needs to know three things:
  //
  // 1. How much weight to assign to each cell (how many particles are in
  // there)
  // 2. How to pack the particles before shipping data around
  // 3. How to unpack the particles after repartitioning
  //
  // Attach the correct functions to the signals inside
  // parallel::distributed::Triangulation, which will be called every time the
  // repartition() or refinement functions are called.
  // These connections only need to be created once, so we might as well
  // have set them up in the constructor of this class, but for the purpose
  // of this example we want to group the particle related instructions.
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
  triangulation.signals.weight.connect(
    [](const typename Triangulation<dim>::cell_iterator &,
       const typename Triangulation<dim>::CellStatus) -> unsigned int {
      return 1000;
    });

  triangulation.signals.weight.connect(
    [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int { return this->cell_weight(cell, status); });

#else
  triangulation.signals.weight.connect(
    [](const typename Triangulation<dim>::cell_iterator &,
       const CellStatus) -> unsigned int { return 1000; });

  triangulation.signals.weight.connect(
    [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator
                        &cell,
        const CellStatus status) -> unsigned int {
      return this->cell_weight(cell, status);
    });
#endif

  if (parameters.model_parameters.sparse_particle_contacts)
    {
      has_sparse_contacts = true;
      sparse_contacts_object.set_parameters(
        parameters.model_parameters.granular_temperature_threshold,
        parameters.model_parameters.solid_fraction_threshold,
        parameters.model_parameters.advect_particles);
    }

  maximum_particle_diameter = 0;
  for (unsigned int particle_type = 0;
       particle_type <
       parameters.lagrangian_physical_properties.particle_type_number;
       particle_type++)
    {
      if (parameters.lagrangian_physical_properties.distribution_type.at(
            particle_type) ==
          Parameters::Lagrangian::SizeDistributionType::uniform)
        {
          size_distribution_object_container[particle_type] =
            std::make_shared<UniformDistribution>(
              parameters.lagrangian_physical_properties
                .particle_average_diameter.at(particle_type));
        }
      else if (parameters.lagrangian_physical_properties.distribution_type.at(
                 particle_type) ==
               Parameters::Lagrangian::SizeDistributionType::normal)
        {
          size_distribution_object_container[particle_type] =
            std::make_shared<NormalDistribution>(
              parameters.lagrangian_physical_properties
                .particle_average_diameter.at(particle_type),
              parameters.lagrangian_physical_properties.particle_size_std.at(
                particle_type),
              parameters.lagrangian_physical_properties
                  .seed_for_distributions[particle_type] +
                this_mpi_process);
        }
      else if (parameters.lagrangian_physical_properties.distribution_type.at(
                 particle_type) ==
               Parameters::Lagrangian::SizeDistributionType::custom)
        {
          size_distribution_object_container[particle_type] =
            std::make_shared<CustomDistribution>(
              parameters.lagrangian_physical_properties.particle_custom_diameter
                .at(particle_type),
              parameters.lagrangian_physical_properties
                .particle_custom_probability.at(particle_type),
              parameters.lagrangian_physical_properties
                  .seed_for_distributions[particle_type] +
                this_mpi_process);
        }
      maximum_particle_diameter = std::max(
        maximum_particle_diameter,
        size_distribution_object_container[particle_type]->find_max_diameter());
    }

  neighborhood_threshold_squared =
    std::pow(parameters.model_parameters.neighborhood_threshold *
               maximum_particle_diameter,
             2);

  if (this_mpi_process == 0)
    input_parameter_inspection(parameters,
                               pcout,
                               size_distribution_object_container);

  grid_motion_object =
    std::make_shared<GridMotion<dim, dim>>(parameters.grid_motion,
                                           simulation_control->get_time_step());

  for (unsigned int i_solid = 0;
       i_solid < parameters.solid_objects->number_solid_surfaces;
       ++i_solid)
    {
      solid_surfaces.push_back(std::make_shared<SerialSolid<dim - 1, dim>>(
        this->parameters.solid_objects->solid_surfaces[i_solid], i_solid));
    }

  for (unsigned int i_solid = 0;
       i_solid < parameters.solid_objects->number_solid_volumes;
       ++i_solid)
    {
      solid_volumes.push_back(std::make_shared<SerialSolid<dim, dim>>(
        this->parameters.solid_objects->solid_volumes[i_solid], i_solid));
    }

  // Generate solid objects
  solid_surfaces_mesh_info.resize(solid_surfaces.size());
  solid_volumes_mesh_info.resize(solid_volumes.size());

  // Resize particle_floating_mesh_in_contact
  if ((solid_surfaces.size() + solid_volumes.size()) > 0)
    {
      has_solid_objects = true;
      contact_manager.particle_floating_mesh_in_contact.resize(
        solid_surfaces.size() + solid_volumes.size());
    }

  // Check if there are periodic boundaries
  for (unsigned int i_bc = 0;
       i_bc < parameters.boundary_conditions.bc_types.size();
       ++i_bc)
    {
      if (parameters.boundary_conditions.bc_types[i_bc] ==
          Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
        {
          has_periodic_boundaries = true;
          break;
        }
    }

  // Assign gravity/acceleration
  g = parameters.lagrangian_physical_properties.g;
}

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
template <int dim>
unsigned int
DEMSolver<dim>::cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const typename parallel::distributed::Triangulation<dim>::CellStatus status)
  const
#else
template <int dim>
unsigned int
DEMSolver<dim>::cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status) const
#endif
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // This determines how important particle work is compared to cell
  // work (by default every cell has a weight of 1000).
  // We set the weight per particle much higher to indicate that
  // the particle load is the only one that is important to distribute
  // in this example. The optimal value of this number depends on the
  // application and can range from 0 (cheap particle operations,
  // expensive cell operations) to much larger than 1000 (expensive
  // particle operations, cheap cell operations, like in this case).
  // This parameter will need to be tuned for the case of DEM.
  const unsigned int particle_weight =
    parameters.model_parameters.load_balance_particle_weight;

  switch (status)
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_PERSIST:
      case parallel::distributed::Triangulation<dim>::CELL_REFINE:

#else
      case CellStatus::cell_will_persist:
      case CellStatus::cell_will_be_refined:
#endif
        // If CELL_PERSIST, do as CELL_REFINE
        {
          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
          break;
        }
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_INVALID:
        break;
#else
      case CellStatus::cell_invalid:
        break;
#endif

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_COARSEN:
#else
      case CellStatus::children_will_be_coarsened:
#endif
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler.n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
          break;
        }

      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
template <int dim>
unsigned int
DEMSolver<dim>::cell_weight_with_mobility_status(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const typename parallel::distributed::Triangulation<dim>::CellStatus status)
  const
#else
template <int dim>
unsigned int
DEMSolver<dim>::cell_weight_with_mobility_status(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status) const
#endif
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  const unsigned int particle_weight =
    parameters.model_parameters.load_balance_particle_weight;

  // Get mobility status of the cell
  const unsigned int cell_mobility_status =
    sparse_contacts_object.check_cell_mobility(cell);

  // Applied a factor on the particle weight regards the mobility status
  // Factor of 1 when mobile cell
  double alpha = 1.0;
  if (cell_mobility_status == sparse_contacts_object.static_active ||
      cell_mobility_status == sparse_contacts_object.advected_active)
    {
      alpha = parameters.model_parameters.active_load_balancing_factor;
    }
  else if (cell_mobility_status == sparse_contacts_object.inactive ||
           cell_mobility_status == sparse_contacts_object.advected)
    {
      alpha = parameters.model_parameters.inactive_load_balancing_factor;
    }

  switch (status)
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_PERSIST:
      case parallel::distributed::Triangulation<dim>::CELL_REFINE:
#else
      case dealii::CellStatus::cell_will_persist:
      case dealii::CellStatus::cell_will_be_refined:

#endif
        {
          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(cell);
          return alpha * n_particles_in_cell * particle_weight;
          break;
        }

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_INVALID:
        break;
#else
      case dealii::CellStatus::cell_invalid:
        break;
#endif

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      case parallel::distributed::Triangulation<dim>::CELL_COARSEN:
#else
      case dealii::CellStatus::children_will_be_coarsened:
#endif
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler.n_particles_in_cell(cell->child(child_index));

          return alpha * n_particles_in_cell * particle_weight;
          break;
        }

      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

template <int dim>
void
DEMSolver<dim>::setup_background_dofs()
{
  FE_Q<dim> background_fe(1);
  background_dh.distribute_dofs(background_fe);

  // Periodic nodes must be mapped otherwise the disabling of contacts will not
  // propagate the mobility status to the periodic nodes during iterations.
  // Identification of periodic nodes is done with the background constraints.
  // Those constraints are not used for any matrix assembly in DEM solver, this
  // approach comes from CFD-DEM coupling where void fraction constraints are
  // used to achieve the periodic node mapping.
  if (has_sparse_contacts && has_periodic_boundaries)
    {
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(background_dh,
                                              locally_relevant_dofs);

      background_constraints.clear();
      background_constraints.reinit(locally_relevant_dofs);

      DoFTools::make_periodicity_constraints(
        background_dh,
        parameters.boundary_conditions.periodic_boundary_0,
        parameters.boundary_conditions.periodic_boundary_1,
        parameters.boundary_conditions.periodic_direction,
        background_constraints);

      background_constraints.close();

      sparse_contacts_object.map_periodic_nodes(background_constraints);
    }
}


template <int dim>
void
DEMSolver<dim>::load_balance()
{
  load_balance_step = load_balance_iteration_check_function();

  if (!load_balance_step)
    return;

  TimerOutput::Scope t(this->computing_timer, "Load balancing");
  // Prepare particle handler for the adaptation of the triangulation to the
  // load
  particle_handler.prepare_for_coarsening_and_refinement();

  pcout << "-->Repartitionning triangulation" << std::endl;
  triangulation.repartition();

  // Unpack the particle handler after the mesh has been repartitioned
  particle_handler.unpack_after_coarsening_and_refinement();

  // Update the container with all the combinations of background and
  // solid cells
  for (unsigned int i_solid = 0; i_solid < solid_surfaces.size(); ++i_solid)
    {
      solid_surfaces_mesh_info[i_solid] =
        solid_surfaces[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  for (unsigned int i_solid = 0; i_solid < solid_volumes.size(); ++i_solid)
    {
      solid_volumes_mesh_info[i_solid] =
        solid_volumes[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  if (has_periodic_boundaries)
    {
      periodic_boundaries_object.map_periodic_cells(
        triangulation, periodic_boundaries_cells_information);

      periodic_offset =
        periodic_boundaries_object.get_periodic_offset_distance();
    }

  // Update neighbors of cells after load balance
  contact_manager.update_cell_neighbors(triangulation,
                                        periodic_boundaries_cells_information,
                                        has_periodic_boundaries,
                                        has_solid_objects);

  boundary_cell_object.build(
    triangulation,
    parameters.floating_walls,
    parameters.boundary_conditions.outlet_boundaries,
    parameters.mesh.check_for_diamond_cells,
    parameters.mesh.expand_particle_wall_contact_search,
    pcout);

  if (parameters.grid_motion.motion_type !=
      Parameters::Lagrangian::GridMotion<dim>::MotionType::none)
    boundary_cell_object.update_boundary_info_after_grid_motion(
      updated_boundary_points_and_normal_vectors);

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
template <int dim>
inline std::function<bool()>
DEMSolver<dim>::set_load_balance_iteration_check_function()
{
  if (parameters.model_parameters.load_balance_method ==
      Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::none)
    {
      return [&] { return false; };
    }
  else if (parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::once)
    {
      return [&] { return check_load_balance_once(); };
    }
  else if (parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::frequent)
    {
      return [&] { return check_load_balance_frequent(); };
    }
  else if (parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::dynamic)
    {
      return [&] { return check_load_balance_dynamic(); };
    }
  else if (parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::
             dynamic_with_sparse_contacts)
    {
      return [&] { return check_load_balance_with_sparse_contacts(); };
    }
  else
    {
      throw std::runtime_error("Specified load balance method is not valid");
    }
}

template <int dim>
inline bool
DEMSolver<dim>::check_load_balance_once()
{
  if (simulation_control->get_step_number() ==
        parameters.model_parameters.load_balance_step ||
      checkpoint_step)
    {
      return true;
    }

  return false;
}

template <int dim>
inline bool
DEMSolver<dim>::check_load_balance_frequent()
{
  if ((simulation_control->get_step_number() %
         parameters.model_parameters.load_balance_frequency ||
       checkpoint_step) == 0)
    {
      return true;
    }

  return false;
}

template <int dim>
inline bool
DEMSolver<dim>::check_load_balance_with_sparse_contacts()
{
  if (simulation_control->get_step_number() %
        parameters.model_parameters.dynamic_load_balance_check_frequency ==
      0)
    {
      // Process to accumulate the load of each process regards the number
      // of cells and particles with their selected weight and with a factor
      // related to the mobility status of the cells
      vector<double> process_to_load_weight(n_mpi_processes, 0.0);

      // Get the particle weight
      const unsigned int particle_weight =
        parameters.model_parameters.load_balance_particle_weight;

      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              // Apply a weight of 1000 to the cell (default value)
              process_to_load_weight[this_mpi_process] += 1000;

              // Get the mobility status of the cell & the number of particles
              const unsigned int cell_mobility_status =
                sparse_contacts_object.check_cell_mobility(cell);
              const unsigned int n_particles_in_cell =
                particle_handler.n_particles_in_cell(cell);

              // Apply a factor on the particle weight regards the
              // mobility status. alpha = 1 by default for mobile cell, but
              // is modified if cell is active or inactive
              double alpha = 1.0;
              if (cell_mobility_status ==
                    sparse_contacts_object.static_active ||
                  cell_mobility_status ==
                    sparse_contacts_object.advected_active)
                {
                  alpha =
                    parameters.model_parameters.active_load_balancing_factor;
                }
              else if (cell_mobility_status ==
                         sparse_contacts_object.inactive ||
                       cell_mobility_status == sparse_contacts_object.advected)
                {
                  alpha =
                    parameters.model_parameters.inactive_load_balancing_factor;
                }

              // Add the particle weight time the number of particles in the
              // cell to the processor load
              process_to_load_weight[this_mpi_process] +=
                alpha * n_particles_in_cell * particle_weight;
            }
        }

      // Exchange information
      double maximum_load_on_proc = 0.0;
      double minimum_load_on_proc = 0.0;
      double total_load           = 0.0;

      maximum_load_on_proc =
        Utilities::MPI::max(*std::max_element(process_to_load_weight.begin(),
                                              process_to_load_weight.end()),
                            mpi_communicator);

      // Find the minimum load on a process
      // First it finds the minimum load on a process, but since values in
      // the vector that are not on this process are 0.0, it looks for
      // values > 1e-8. After that, it finds the minimum load of all the
      // processors
      minimum_load_on_proc =
        Utilities::MPI::min(*std::min_element(process_to_load_weight.begin(),
                                              process_to_load_weight.end(),
                                              [](double a, double b) {
                                                return (a > 1e-8) ?
                                                         (b > 1e-8 ? a < b :
                                                                     true) :
                                                         false;
                                              }),
                            mpi_communicator);

      // Get the total load
      total_load =
        Utilities::MPI::sum(std::accumulate(process_to_load_weight.begin(),
                                            process_to_load_weight.end(),
                                            0.0),
                            mpi_communicator);

      if ((maximum_load_on_proc - minimum_load_on_proc) >
            parameters.model_parameters.load_balance_threshold *
              (total_load / n_mpi_processes) ||
          checkpoint_step)
        {
          return true;
        }
    }

  // Clear and connect a new cell weight function
  triangulation.signals.weight.disconnect_all_slots();

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
  triangulation.signals.weight.connect(
    [](const typename Triangulation<dim>::cell_iterator &,
       const typename Triangulation<dim>::CellStatus) -> unsigned int {
      return 1000;
    });

  triangulation.signals.weight.connect(
    [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int {
      return this->cell_weight_with_mobility_status(cell, status);
    });

#else
  triangulation.signals.weight.connect(
    [](const typename Triangulation<dim>::cell_iterator &,
       const CellStatus) -> unsigned int { return 1000; });

  triangulation.signals.weight.connect(
    [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator
                        &cell,
        const CellStatus status) -> unsigned int {
      return this->cell_weight_with_mobility_status(cell, status);
    });
#endif

  return false;
}

template <int dim>
inline bool
DEMSolver<dim>::check_load_balance_dynamic()
{
  if (simulation_control->get_step_number() %
        parameters.model_parameters.dynamic_load_balance_check_frequency ==
      0)
    {
      unsigned int maximum_particle_number_on_proc = 0;
      unsigned int minimum_particle_number_on_proc = 0;

      maximum_particle_number_on_proc =
        Utilities::MPI::max(particle_handler.n_locally_owned_particles(),
                            mpi_communicator);
      minimum_particle_number_on_proc =
        Utilities::MPI::min(particle_handler.n_locally_owned_particles(),
                            mpi_communicator);

      // Execute load balancing if difference of load between processors is
      // larger than threshold of the load per processor
      if ((maximum_particle_number_on_proc - minimum_particle_number_on_proc) >
            parameters.model_parameters.load_balance_threshold *
              (particle_handler.n_global_particles() / n_mpi_processes) ||
          checkpoint_step)
        {
          return true;
        }
    }

  return false;
}

template <int dim>
inline std::function<bool()>
DEMSolver<dim>::set_contact_search_iteration_function()
{
  if (parameters.model_parameters.contact_detection_method ==
      Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::constant)
    {
      return [&] { return check_contact_search_iteration_constant(); };
    }
  else if (parameters.model_parameters.contact_detection_method ==
           Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::
             dynamic)
    {
      return [&] { return check_contact_search_iteration_dynamic(); };
    }
  else
    {
      throw std::runtime_error(
        "Specified contact detection method is not valid");
    }
}

template <int dim>
inline bool
DEMSolver<dim>::check_contact_search_iteration_dynamic()
{
  bool sorting_in_subdomains_step =
    (particles_insertion_step || load_balance_step || contact_detection_step);

  return find_particle_contact_detection_step<dim>(
    particle_handler,
    simulation_control->get_time_step(),
    smallest_contact_search_criterion,
    mpi_communicator,
    sorting_in_subdomains_step,
    displacement,
    (simulation_control->get_step_number() % contact_detection_frequency) == 0);
}

template <int dim>
inline bool
DEMSolver<dim>::check_contact_search_iteration_constant()
{
  return (
    (simulation_control->get_step_number() % contact_detection_frequency) == 0);
}

template <int dim>
bool
DEMSolver<dim>::insert_particles()
{
  if ((simulation_control->get_step_number() % insertion_frequency) == 1)
    {
      insertion_object->insert(particle_handler, triangulation, parameters);
      return true;
    }
  return false;
}

template <int dim>
void
DEMSolver<dim>::update_moment_of_inertia(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  std::vector<double>                     &MOI)
{
  MOI.resize(torque.size());

  for (auto &particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      MOI[particle.get_local_index()] =
        0.1 * particle_properties[DEM::PropertiesIndex::mass] *
        particle_properties[DEM::PropertiesIndex::dp] *
        particle_properties[DEM::PropertiesIndex::dp];
    }
}

template <int dim>
void
DEMSolver<dim>::particle_wall_contact_force()
{
  // Particle-wall contact force
  particle_wall_contact_force_object->calculate_particle_wall_contact_force(
    contact_manager.particle_wall_in_contact,
    simulation_control->get_time_step(),
    torque,
    force);

  if (parameters.forces_torques.calculate_force_torque)
    {
      forces_boundary_information[simulation_control->get_step_number()] =
        particle_wall_contact_force_object->get_force();
      torques_boundary_information[simulation_control->get_step_number()] =
        particle_wall_contact_force_object->get_torque();
    }

  // Particle-floating wall contact force
  if (parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_contact_force_object->calculate_particle_wall_contact_force(
        contact_manager.particle_floating_wall_in_contact,
        simulation_control->get_time_step(),
        torque,
        force);
    }

  // Particle - floating mesh contact force
  if (has_solid_objects)
    {
      particle_wall_contact_force_object
        ->calculate_particle_floating_wall_contact_force(
          contact_manager.particle_floating_mesh_in_contact,
          simulation_control->get_time_step(),
          torque,
          force,
          solid_surfaces);
    }

  particle_point_line_contact_force_object
    .calculate_particle_point_contact_force(
      &contact_manager.particle_points_in_contact,
      parameters.lagrangian_physical_properties,
      force);

  if (dim == 3)
    {
      particle_point_line_contact_force_object
        .calculate_particle_line_contact_force(
          &contact_manager.particle_lines_in_contact,
          parameters.lagrangian_physical_properties,
          force);
    }
}

template <int dim>
void
DEMSolver<dim>::finish_simulation()
{
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    {
      this->computing_timer.print_summary();
    }

  // Testing
  if (parameters.test.enabled)
    {
      switch (parameters.test.test_type)
        {
          case Parameters::Testing::TestType::particles:
            {
              visualization_object.print_xyz(particle_handler,
                                             mpi_communicator,
                                             pcout);
              break;
            }
          case Parameters::Testing::TestType::mobility_status:
            {
              // Get mobility status vector sorted by cell id
              Vector<float> mobility_status(triangulation.n_active_cells());
              sparse_contacts_object.get_mobility_status_vector(
                mobility_status);

              // Output mobility status vector
              visualization_object.print_intermediate_format(mobility_status,
                                                             background_dh,
                                                             mpi_communicator);
              break;
            }
          case Parameters::Testing::TestType::subdomain:
            {
              // Get mobility status vector sorted by cell id
              Vector<float> subdomain(triangulation.n_active_cells());
              for (unsigned int i = 0; i < subdomain.size(); ++i)
                subdomain(i) = triangulation.locally_owned_subdomain();

              // Output subdomain vector
              visualization_object.print_intermediate_format(subdomain,
                                                             background_dh,
                                                             mpi_communicator);
              break;
            }
        }
    }

  // Outputting force and torques over boundary
  if (parameters.forces_torques.calculate_force_torque)
    {
      write_forces_torques_output_results<dim>(
        parameters.forces_torques.force_torque_output_name,
        parameters.forces_torques.output_frequency,
        triangulation.get_boundary_ids(),
        simulation_control->get_time_step(),
        forces_boundary_information,
        torques_boundary_information);
    }
}

template <int dim>
std::shared_ptr<Insertion<dim>>
DEMSolver<dim>::set_insertion_type(const DEMSolverParameters<dim> &parameters)
{
  if (parameters.insertion_info.insertion_method ==
      Parameters::Lagrangian::InsertionInfo::InsertionMethod::file)
    {
      insertion_object =
        std::make_shared<InsertionFile<dim>>(size_distribution_object_container,
                                             triangulation,
                                             parameters);
    }
  else if (parameters.insertion_info.insertion_method ==
           Parameters::Lagrangian::InsertionInfo::InsertionMethod::list)
    {
      insertion_object =
        std::make_shared<InsertionList<dim>>(size_distribution_object_container,
                                             triangulation,
                                             parameters);
    }
  else if (parameters.insertion_info.insertion_method ==
           Parameters::Lagrangian::InsertionInfo::InsertionMethod::plane)
    {
      insertion_object = std::make_shared<InsertionPlane<dim>>(
        size_distribution_object_container, triangulation, parameters);
    }
  else if (parameters.insertion_info.insertion_method ==
           Parameters::Lagrangian::InsertionInfo::InsertionMethod::volume)
    {
      insertion_object = std::make_shared<InsertionVolume<dim>>(
        size_distribution_object_container,
        triangulation,
        parameters,
        maximum_particle_diameter);
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
  if (parameters.model_parameters.integration_method ==
      Parameters::Lagrangian::ModelParameters::IntegrationMethod::
        velocity_verlet)
    {
      integrator_object = std::make_shared<VelocityVerletIntegrator<dim>>();
    }
  else if (parameters.model_parameters.integration_method ==
           Parameters::Lagrangian::ModelParameters::IntegrationMethod::
             explicit_euler)
    {
      integrator_object = std::make_shared<ExplicitEulerIntegrator<dim>>();
    }
  else if (parameters.model_parameters.integration_method ==
           Parameters::Lagrangian::ModelParameters::IntegrationMethod::gear3)
    {
      integrator_object = std::make_shared<Gear3Integrator<dim>>();
    }
  else
    {
      throw "The chosen integration method is invalid";
    }
  return integrator_object;
}

template <int dim>
void
DEMSolver<dim>::write_output_results()
{
  TimerOutput::Scope t(this->computing_timer, "Output VTU");

  const std::string folder = parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    parameters.simulation_control.output_name;
  const unsigned int iter        = simulation_control->get_step_number();
  const double       time        = simulation_control->get_current_time();
  const unsigned int group_files = parameters.simulation_control.group_files;

  // Write particles
  Visualization<dim> particle_data_out;
  particle_data_out.build_patches(particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler,
                            particle_data_out,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            mpi_communicator);

  if (simulation_control->get_output_boundaries())
    {
      DataOutFaces<dim> data_out_faces;

      // Setup background dofs to initiate right sized boundary_id vector
      Vector<float> boundary_id(background_dh.n_dofs());

      // Attach the boundary_id to data_out_faces object
      BoundaryPostprocessor<dim> boundary;
      data_out_faces.attach_dof_handler(background_dh);
      data_out_faces.add_data_vector(boundary_id, boundary);
      data_out_faces.build_patches();

      write_boundaries_vtu<dim>(
        data_out_faces, folder, time, iter, this->mpi_communicator);
    }
  if (parameters.post_processing.force_chains)
    {
      // Force chains visualization
      particles_force_chains_object =
        set_force_chains_contact_force_model(parameters);
      particles_force_chains_object->calculate_force_chains(contact_manager);
      particles_force_chains_object->write_force_chains(
        parameters,
        particles_pvdhandler_force_chains,
        this->mpi_communicator,
        folder,
        iter,
        time);
    }

  // Write all solid objects
  for (const auto &solid_object : solid_surfaces)
    solid_object->write_output_results(simulation_control);

  for (const auto &solid_object : solid_volumes)
    solid_object->write_output_results(simulation_control);
}

template <int dim>
void
DEMSolver<dim>::post_process_results()
{
  post_processing_object.write_post_processing_results(
    triangulation,
    grid_pvdhandler,
    background_dh,
    particle_handler,
    parameters,
    simulation_control->get_current_time(),
    simulation_control->get_step_number(),
    mpi_communicator,
    sparse_contacts_object);
}

template <int dim>
void
DEMSolver<dim>::report_statistics()
{
  // Update statistics on contact list
  double number_of_list_built_since_last_log =
    double(contact_build_number) - contact_list.total;
  contact_list.max =
    std::max(number_of_list_built_since_last_log, contact_list.max);
  contact_list.min =
    std::min(number_of_list_built_since_last_log, contact_list.min);
  contact_list.total += number_of_list_built_since_last_log;
  contact_list.average = contact_list.total /
                         (simulation_control->get_step_number()) *
                         simulation_control->get_log_frequency();

  // Calculate statistics on the particles
  statistics translational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::dem_statistic_variable::translational_kinetic_energy>(
    particle_handler, mpi_communicator);
  statistics rotational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::dem_statistic_variable::rotational_kinetic_energy>(particle_handler,
                                                            mpi_communicator);
  statistics velocity =
    calculate_granular_statistics<dim, DEM::dem_statistic_variable::velocity>(
      particle_handler, mpi_communicator);
  statistics omega =
    calculate_granular_statistics<dim, DEM::dem_statistic_variable::omega>(
      particle_handler, mpi_communicator);

  if (this_mpi_process == 0)
    {
      TableHandler report;

      report.declare_column("Variable");
      report.declare_column("Min");
      report.declare_column("Max");
      report.declare_column("Average");
      report.declare_column("Total");
      add_statistics_to_table_handler("Contact list generation",
                                      contact_list,
                                      report);
      add_statistics_to_table_handler("Velocity magnitude", velocity, report);
      add_statistics_to_table_handler("Angular velocity magnitude",
                                      omega,
                                      report);
      add_statistics_to_table_handler("Translational kinetic energy",
                                      translational_kinetic_energy,
                                      report);
      add_statistics_to_table_handler("Rotational kinetic energy",
                                      rotational_kinetic_energy,
                                      report);



      report.set_scientific("Min", true);
      report.set_scientific("Max", true);
      report.set_scientific("Average", true);
      report.set_scientific("Total", true);

      report.write_text(std::cout, dealii::TableHandler::org_mode_table);
    }

  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::iteration)
    {
      this->computing_timer.print_summary();
    }
}

template <int dim>
void
DEMSolver<dim>::solve()
{
  // Reading mesh
  read_mesh(parameters.mesh,
            parameters.restart.restart,
            pcout,
            triangulation,
            parameters.boundary_conditions);

  // Store information about floating mesh/background mesh intersection
  for (unsigned int i_solid = 0; i_solid < solid_surfaces.size(); ++i_solid)
    {
      solid_surfaces_mesh_info[i_solid] =
        solid_surfaces[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  for (unsigned int i_solid = 0; i_solid < solid_volumes.size(); ++i_solid)
    {
      solid_volumes_mesh_info[i_solid] =
        solid_volumes[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  // Set insertion object type before the restart because the restart only
  // rebuilds the member of the insertion object
  insertion_object = set_insertion_type(parameters);

  load_balance_iteration_check_function =
    set_load_balance_iteration_check_function();

  contact_detection_iteration_check_function =
    set_contact_search_iteration_function();

  if (parameters.restart.restart == true)
    {
      read_checkpoint(computing_timer,
                      parameters,
                      simulation_control,
                      particles_pvdhandler,
                      grid_pvdhandler,
                      triangulation,
                      particle_handler,
                      insertion_object,
                      solid_surfaces);

      displacement.resize(particle_handler.get_max_local_particle_index());
      force.resize(displacement.size());
      torque.resize(displacement.size());

      update_moment_of_inertia(particle_handler, MOI);

      checkpoint_step = true;
    }

  // Store information about floating mesh/background mesh intersection
  for (unsigned int i_solid = 0; i_solid < solid_surfaces.size(); ++i_solid)
    {
      solid_surfaces_mesh_info[i_solid] =
        solid_surfaces[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  for (unsigned int i_solid = 0; i_solid < solid_volumes.size(); ++i_solid)
    {
      solid_volumes_mesh_info[i_solid] =
        solid_volumes[i_solid]->map_solid_in_background_triangulation(
          triangulation);
    }

  // Find the smallest contact search frequency criterion between (smallest
  // cell size - largest particle radius) and (security factor * (blob
  // diameter - 1) *  largest particle radius). This value is used in
  // find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(triangulation) -
              maximum_particle_diameter * 0.5),
             (parameters.model_parameters.dynamic_contact_search_factor *
              (parameters.model_parameters.neighborhood_threshold - 1) *
              maximum_particle_diameter * 0.5));

  // Find the smallest cell size and use this as the floating mesh mapping
  // criterion

  double mapping_criterion_constant;
  if constexpr (dim == 2)
    mapping_criterion_constant = 0.70710678118; // 2^-0.5

  if constexpr (dim == 3)
    mapping_criterion_constant = 0.57735026919; // 3^-0.5

  smallest_floating_mesh_mapping_criterion =
    mapping_criterion_constant *
    GridTools::minimal_cell_diameter(triangulation);

  // The edge case comes when the cell are completely square/cubic. In that
  // case, every sides of a cell are 2^-0.5 or 3^-0.5 times the cell_diameter.
  // We want to refresh the mapping each time the solid-objet pass through a
  // cell or there will be late contact detection. Thus, we use this value.


  if (has_periodic_boundaries)
    {
      periodic_boundaries_object.set_periodic_boundaries_information(
        parameters.boundary_conditions.periodic_boundary_0,
        parameters.boundary_conditions.periodic_boundary_1,
        parameters.boundary_conditions.periodic_direction);

      periodic_boundaries_object.map_periodic_cells(
        triangulation, periodic_boundaries_cells_information);

      // Temporary offset calculation : works only for one set of periodic
      // boundary on an axis.
      periodic_offset =
        periodic_boundaries_object.get_periodic_offset_distance();
    }

  // Find cell neighbors
  contact_manager.execute_cell_neighbors_search(
    triangulation,
    periodic_boundaries_cells_information,
    has_periodic_boundaries,
    has_solid_objects);

  // Finding boundary cells with faces
  boundary_cell_object.build(
    triangulation,
    parameters.floating_walls,
    parameters.boundary_conditions.outlet_boundaries,
    parameters.mesh.check_for_diamond_cells,
    parameters.mesh.expand_particle_wall_contact_search,
    pcout);

  // Setting chosen contact force, insertion and integration methods
  integrator_object = set_integrator_type(parameters);
  particle_particle_contact_force_object =
    set_particle_particle_contact_force_model(parameters);
  particle_wall_contact_force_object =
    set_particle_wall_contact_force_model(parameters, triangulation);

  // Setup background dof
  setup_background_dofs();

  // DEM engine iterator:
  while (simulation_control->integrate())
    {
      simulation_control->print_progression(pcout);
      if (simulation_control->is_verbose_iteration())
        report_statistics();

      // Grid motion
      if (parameters.grid_motion.motion_type !=
          Parameters::Lagrangian::GridMotion<dim>::MotionType::none)
        {
          grid_motion_object->move_grid(triangulation);
          boundary_cell_object.update_boundary_info_after_grid_motion(
            updated_boundary_points_and_normal_vectors);
        }

      // Keep track if particles were inserted this step
      particles_insertion_step = insert_particles();

      if (particles_insertion_step)
        displacement.resize(particle_handler.get_max_local_particle_index());

      // Load balancing
      load_balance();

      if (load_balance_step || checkpoint_step)
        {
          displacement.resize(particle_handler.get_max_local_particle_index());

          if (has_sparse_contacts)
            {
              sparse_contacts_object.update_local_and_ghost_cell_set(
                background_dh);
            }
        }

      // Check to see if it is contact search step
      contact_detection_step = contact_detection_iteration_check_function();

      bool solid_object_map_step = false;
      // Check to see if floating meshes need to be mapped in background mesh
      if (has_solid_objects)
        {
          solid_object_map_step = find_floating_mesh_mapping_step(
            smallest_floating_mesh_mapping_criterion, this->solid_surfaces);

          if (solid_object_map_step)
            {
              // Update floating mesh information in the container manager
              for (unsigned int i_solid = 0; i_solid < solid_surfaces.size();
                   ++i_solid)
                {
                  solid_surfaces_mesh_info[i_solid] =
                    solid_surfaces[i_solid]
                      ->map_solid_in_background_triangulation(triangulation);
                }

              for (unsigned int i_solid = 0; i_solid < solid_volumes.size();
                   ++i_solid)
                {
                }
            }
        }

      contact_detection_step = contact_detection_step || solid_object_map_step;

      // Sort particles in cells
      if (particles_insertion_step || load_balance_step ||
          contact_detection_step || checkpoint_step)
        {
          // Particles displacement if passing through a periodic boundary
          periodic_boundaries_object.execute_particles_displacement(
            particle_handler, periodic_boundaries_cells_information);

          particle_handler.sort_particles_into_subdomains_and_cells();

          if (has_sparse_contacts && !simulation_control->is_at_start())
            {
              // Compute cell mobility for all cells after the sort particles
              // into subdomains and cells and exchange ghost particles
              sparse_contacts_object.identify_mobility_status(
                background_dh,
                particle_handler,
                triangulation.n_active_cells(),
                mpi_communicator);
            }

          displacement.resize(particle_handler.get_max_local_particle_index());
          force.resize(displacement.size());
          torque.resize(displacement.size());

          // Updating moment of inertia container
          update_moment_of_inertia(particle_handler, MOI);

          particle_handler.exchange_ghost_particles(true);

          // Reset checkpoint step
          checkpoint_step = false;

          // Execute broad search by filling containers of particle-particle
          // contact pair candidates and containers of particle-wall
          // contact pair candidates
          if (!(has_sparse_contacts && contact_build_number > 1))
            {
              contact_manager.execute_particle_particle_broad_search(
                particle_handler, has_periodic_boundaries);

              contact_manager.execute_particle_wall_broad_search(
                particle_handler,
                boundary_cell_object,
                solid_surfaces_mesh_info,
                parameters.floating_walls,
                simulation_control->get_current_time(),
                has_solid_objects);
            }
          else
            {
              contact_manager.execute_particle_particle_broad_search(
                particle_handler,
                sparse_contacts_object,
                has_periodic_boundaries);

              contact_manager.execute_particle_wall_broad_search(
                particle_handler,
                boundary_cell_object,
                solid_surfaces_mesh_info,
                parameters.floating_walls,
                simulation_control->get_current_time(),
                sparse_contacts_object,
                has_solid_objects);
            }

          // Updating number of contact builds
          contact_build_number++;

          // Update contacts, remove replicates and add new contact pairs
          // to the contact containers when particles are exchanged between
          // processors
          contact_manager.update_contacts(has_periodic_boundaries);

          // Updates the iterators to particles in local-local contact
          // containers
          contact_manager.update_local_particles_in_cells(
            particle_handler, load_balance_step, has_periodic_boundaries);

          // Execute fine search by updating particle-particle contact
          // containers according to the neighborhood threshold
          contact_manager.execute_particle_particle_fine_search(
            neighborhood_threshold_squared,
            has_periodic_boundaries,
            periodic_offset);

          // Execute fine search by updating particle-wall contact containers
          // according to the neighborhood threshold
          contact_manager.execute_particle_wall_fine_search(
            parameters.floating_walls,
            simulation_control->get_current_time(),
            neighborhood_threshold_squared,
            has_solid_objects);
        }
      else
        {
          particle_handler.update_ghost_particles();
        }

      // Particle-particle contact force
      particle_particle_contact_force_object
        ->calculate_particle_particle_contact_force(
          contact_manager,
          simulation_control->get_time_step(),
          torque,
          force,
          periodic_offset);

      // We have to update the positions of the points on boundary faces and
      // their normal vectors here. The update_contacts deletes the
      // particle-wall contact candidate if it exists in the contact list. As
      // a result, when we update the points on boundary faces and their
      // normal vectors, update_contacts deletes it from the output of broad
      // search and they are not updated in the contact force calculations.
      if (parameters.grid_motion.motion_type !=
          Parameters::Lagrangian::GridMotion<dim>::MotionType::none)
        grid_motion_object
          ->update_boundary_points_and_normal_vectors_in_contact_list(
            contact_manager.particle_wall_in_contact,
            updated_boundary_points_and_normal_vectors);

      // Move the solid triangulations, previous time must be used here
      // instead of current time.
      for (auto &solid_object : solid_surfaces)
        solid_object->move_solid_triangulation(
          simulation_control->get_time_step(),
          simulation_control->get_previous_time());

      for (auto &solid_object : solid_volumes)
        solid_object->move_solid_triangulation(
          simulation_control->get_time_step(),
          simulation_control->get_previous_time());

      // Particles-walls contact force:
      particle_wall_contact_force();

      // Integration correction step (after force calculation)
      // In the first step, we have to obtain location of particles at
      // half-step time
      if (simulation_control->get_step_number() == 0)
        {
          integrator_object->integrate_half_step_location(
            particle_handler,
            g,
            simulation_control->get_time_step(),
            torque,
            force,
            MOI);
        }
      else // Step number > 0
        {
          if (!(has_sparse_contacts && contact_build_number > 1))
            {
              integrator_object->integrate(particle_handler,
                                           g,
                                           simulation_control->get_time_step(),
                                           torque,
                                           force,
                                           MOI);
            }
          else
            {
              integrator_object->integrate(particle_handler,
                                           g,
                                           simulation_control->get_time_step(),
                                           torque,
                                           force,
                                           MOI,
                                           triangulation,
                                           sparse_contacts_object);
            }
        }

      // Visualization
      if (simulation_control->is_output_iteration())
        {
          write_output_results();
        }

      if (parameters.forces_torques.calculate_force_torque &&
          (this_mpi_process == 0) &&
          (simulation_control->get_step_number() %
             parameters.forces_torques.output_frequency ==
           0) &&
          (parameters.forces_torques.force_torque_verbosity ==
           Parameters::Verbosity::verbose))
        write_forces_torques_output_locally(
          forces_boundary_information[simulation_control->get_step_number()],
          torques_boundary_information[simulation_control->get_step_number()]);

      // Post-processing
      if (parameters.post_processing.Lagrangian_post_processing &&
          simulation_control->is_output_iteration())
        {
          post_process_results();
        }

      if (parameters.restart.checkpoint &&
          simulation_control->get_step_number() %
              parameters.restart.frequency ==
            0)
        {
          write_checkpoint(computing_timer,
                           parameters,
                           simulation_control,
                           particles_pvdhandler,
                           grid_pvdhandler,
                           triangulation,
                           particle_handler,
                           insertion_object,
                           solid_surfaces,
                           pcout,
                           mpi_communicator);
        }
    }

  finish_simulation();
}

template class DEMSolver<2>;
template class DEMSolver<3>;
