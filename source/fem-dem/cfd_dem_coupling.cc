#include <core/grids.h>
#include <core/solutions_output.h>

#include <dem/explicit_euler_integrator.h>
#include <dem/gear3_integrator.h>
#include <dem/post_processing.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <fem-dem/cfd_dem_coupling.h>

#include <fstream>
#include <sstream>

template <int dim>
bool
check_contact_detection_method(
  unsigned int                          counter,
  CFDDEMSimulationParameters<dim>      &param,
  std::vector<double>                  &displacement,
  Particles::ParticleHandler<dim, dim> &particle_handler,
  MPI_Comm                              mpi_communicator,
  std::shared_ptr<SimulationControl>    simulation_control,
  bool                                  contact_detection_step,
  bool                                  checkpoint_step,
  bool                                  load_balance_step,
  double                                smallest_contact_search_criterion)
{
  if (param.dem_parameters.model_parameters.contact_detection_method ==
      Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::constant)
    {
      return ((counter % param.dem_parameters.model_parameters
                           .contact_detection_frequency) == 0);
    }
  else if (param.dem_parameters.model_parameters.contact_detection_method ==
           Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::
             dynamic)
    {
      // The sorting into subdomain step checks whether or not the current time
      // step is a step that requires sorting particles into subdomains and
      // cells.
      // This is applicable if any of the following three conditions apply:if
      // its a load balancing step, a restart simulation step, or a contact
      // detection tsep.
      bool sorting_in_subdomains_step =
        (checkpoint_step || load_balance_step || contact_detection_step);

      if (sorting_in_subdomains_step)
        displacement.resize(particle_handler.get_max_local_particle_index());

      contact_detection_step = find_particle_contact_detection_step<dim>(
        particle_handler,
        simulation_control->get_time_step() / param.cfd_dem.coupling_frequency,
        smallest_contact_search_criterion,
        mpi_communicator,
        sorting_in_subdomains_step,
        displacement);

      return contact_detection_step;
    }
  else
    {
      throw std::runtime_error(
        "Specified contact detection method is not valid");
    }
}

template <int dim>
bool
check_load_balance_method(
  CFDDEMSimulationParameters<dim>      &param,
  Particles::ParticleHandler<dim, dim> &particle_handler,
  const MPI_Comm                       &mpi_communicator,
  const unsigned int                    n_mpi_processes,
  std::shared_ptr<SimulationControl>    simulation_control)
{ // Setting load-balance method (single-step, frequent or dynamic)
  if (param.dem_parameters.model_parameters.load_balance_method ==
      Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::once)
    {
      return (simulation_control->get_step_number() ==
              param.dem_parameters.model_parameters.load_balance_step);
    }
  else if (param.dem_parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::frequent)
    {
      return (simulation_control->get_step_number() %
                param.dem_parameters.model_parameters.load_balance_frequency ==
              0);
    }
  else if (param.dem_parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::dynamic)
    {
      bool load_balance_step = false;
      if (simulation_control->get_step_number() %
            param.dem_parameters.model_parameters
              .dynamic_load_balance_check_frequency ==
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

          if ((maximum_particle_number_on_proc -
               minimum_particle_number_on_proc) >
              param.dem_parameters.model_parameters.load_balance_threshold *
                (particle_handler.n_global_particles() / n_mpi_processes))
            {
              load_balance_step = true;
            }
        }
      return load_balance_step;
    }
  else if (param.dem_parameters.model_parameters.load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::none)
    {
      return false;
    }
  else
    {
      throw std::runtime_error("Specified load balance method is not valid");
    }
}

// Constructor for class CFD-DEM
template <int dim>
CFDDEMSolver<dim>::CFDDEMSolver(CFDDEMSimulationParameters<dim> &nsparam)
  : GLSVANSSolver<dim>(nsparam)
  , has_periodic_boundaries(false)
  , has_sparse_contacts(false)
  , this_mpi_process(Utilities::MPI::this_mpi_process(this->mpi_communicator))
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(this->mpi_communicator))
{}

template <int dim>
CFDDEMSolver<dim>::~CFDDEMSolver()
{}

template <int dim>
void
CFDDEMSolver<dim>::read_dem()
{
  this->pcout << "Reading DEM checkpoint " << std::endl;

  std::string prefix =
    this->cfd_dem_simulation_parameters.void_fraction->dem_file_name;

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> this->particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      try
        {
          parallel_triangulation->load(filename.c_str());
        }
      catch (...)
        {
          AssertThrow(false,
                      ExcMessage("Cannot open snapshot mesh file or read the"
                                 "triangulation stored there."));
        }
    }
  else
    {
      throw std::runtime_error(
        "VANS equations currently do not support triangulations other than parallel::distributed");
    }

  // Deserialize particles have the triangulation has been read
  this->particle_handler.deserialize();

  this->pcout << "Finished reading DEM checkpoint " << std::endl
              << this->particle_handler.n_global_particles()
              << " particles are in the simulation " << std::endl;

  write_DEM_output_results();
}

template <int dim>
void
CFDDEMSolver<dim>::write_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "write_checkpoint");

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  std::string prefix_grid      = prefix + "_postprocess_data";
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);
      this->pvdhandler.save(prefix);
      particles_pvdhandler.save(prefix_particles);

      if (this->dem_parameters.post_processing.Lagrangian_post_processing)
        grid_pvdhandler.save(prefix_grid);

      if (this->simulation_parameters.flow_control.enable_flow_control)
        this->flow_control.save(prefix);
    }

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << this->particle_handler;

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;

  std::vector<const GlobalVectorType *> sol_set_transfer;
  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<const GlobalVectorType *> av_set_transfer =
        this->average_velocities->save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    system_trans_vectors(this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  // Prepare particle handler for serialization
  this->particle_handler.prepare_for_serialization();

  // Void Fraction
  std::vector<const GlobalVectorType *> vf_set_transfer;
  vf_set_transfer.push_back(&this->nodal_void_fraction_relevant);
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_set_transfer.push_back(&this->previous_void_fraction[i]);
    }

  this->multiphysics->write_checkpoint();

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    vf_system_trans_vectors(this->void_fraction_dof_handler);
  vf_system_trans_vectors.prepare_for_serialization(vf_set_transfer);

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      std::string triangulationName = prefix + ".triangulation";
      parallel_triangulation->save(prefix + ".triangulation");
    }
}

template <int dim>
void
CFDDEMSolver<dim>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "read_checkpoint");
  std::string        prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  std::string prefix_grid      = prefix + "_postprocess_data";

  this->simulation_control->read(prefix);
  this->pvdhandler.read(prefix);
  particles_pvdhandler.read(prefix_particles);

  if (this->dem_parameters.post_processing.Lagrangian_post_processing)
    {
      grid_pvdhandler.read(prefix_grid);
    }

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> this->particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  try
    {
      if (auto parallel_triangulation =
            dynamic_cast<parallel::distributed::Triangulation<dim> *>(
              this->triangulation.get()))
        parallel_triangulation->load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }

  this->setup_dofs();

  // Remap periodic nodes after setup of dofs
  if (has_periodic_boundaries && has_sparse_contacts)
    {
      sparse_contacts_object.map_periodic_nodes(
        this->void_fraction_constraints);
    }

  // Velocity Vectors
  std::vector<GlobalVectorType *> x_system(1 + this->previous_solutions.size());

  GlobalVectorType distributed_system(this->locally_owned_dofs,
                                      this->mpi_communicator);

  x_system[0] = &(distributed_system);

  std::vector<GlobalVectorType> distributed_previous_solutions;

  distributed_previous_solutions.reserve(this->previous_solutions.size());

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(this->locally_owned_dofs, this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }

  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    system_trans_vectors(this->dof_handler);

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<GlobalVectorType *> sum_vectors =
        this->average_velocities->read(prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);

  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      this->previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Void Fraction Vectors
  std::vector<GlobalVectorType *> vf_system(
    1 + this->previous_void_fraction.size());

  GlobalVectorType vf_distributed_system(this->locally_owned_dofs_voidfraction,
                                         this->mpi_communicator);

  vf_system[0] = &(vf_distributed_system);

  std::vector<GlobalVectorType> vf_distributed_previous_solutions;

  vf_distributed_previous_solutions.reserve(
    this->previous_void_fraction.size());

  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_distributed_previous_solutions.emplace_back(
        GlobalVectorType(this->locally_owned_dofs_voidfraction,
                         this->mpi_communicator));
      vf_system[i + 1] = &vf_distributed_previous_solutions[i];
    }

  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    vf_system_trans_vectors(this->void_fraction_dof_handler);

  vf_system_trans_vectors.deserialize(vf_system);

  this->nodal_void_fraction_relevant = vf_distributed_system;
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      this->previous_void_fraction[i] = vf_distributed_previous_solutions[i];
    }

  if (this->simulation_parameters.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);
    }

  this->multiphysics->read_checkpoint();

  // Deserialize particles have the triangulation has been read
  this->particle_handler.deserialize();
}

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
template <int dim>
unsigned int
CFDDEMSolver<dim>::cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const typename parallel::distributed::Triangulation<dim>::CellStatus status)
  const
#else
template <int dim>
unsigned int
CFDDEMSolver<dim>::cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status) const
#endif


{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // This determines how important particle work is compared to cell
  // work (by default every cell has a weight of 1000).
  // We set the weight per particle higher to indicate that
  // the particle load is more important than the fluid load. The optimal
  // value of this number depends on the application and can range from 0
  // (cheap particle operations, expensive cell operations) to much larger
  // than 1000 (expensive particle operations, cheap cell operations, like in
  // this case). This parameter will need to be tuned for different cases of
  // CFD-DEM coupling.
  const unsigned int particle_weight =
    dem_parameters.model_parameters.load_balance_particle_weight;

  // This does not use adaptive refinement, therefore every cell
  // should have the status CELL_PERSIST. However this function can also
  // be used to distribute load during refinement, therefore we consider
  // refined or coarsened cells as well.

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
  if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
      status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
#else
  if (status == CellStatus::cell_will_persist ||
      status == CellStatus::cell_will_be_refined)
#endif

    {
      const unsigned int n_particles_in_cell =
        this->particle_handler.n_particles_in_cell(cell);
      return n_particles_in_cell * particle_weight;
    }
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
  else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
#else
  else if (status == CellStatus::children_will_be_coarsened)
#endif
    {
      unsigned int n_particles_in_cell = 0;

      for (unsigned int child_index = 0;
           child_index < GeometryInfo<dim>::max_children_per_cell;
           ++child_index)
        n_particles_in_cell +=
          this->particle_handler.n_particles_in_cell(cell->child(child_index));

      return n_particles_in_cell * particle_weight;
    }

  Assert(false, ExcInternalError());
  return 0;
}

template <int dim>
void
CFDDEMSolver<dim>::load_balance()
{
  std::vector<const GlobalVectorType *> sol_set_transfer;
  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    system_trans_vectors(this->dof_handler);
  system_trans_vectors.prepare_for_coarsening_and_refinement(sol_set_transfer);


  // Void Fraction
  std::vector<const GlobalVectorType *> vf_set_transfer;
  vf_set_transfer.push_back(&this->nodal_void_fraction_relevant);
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_set_transfer.push_back(&this->previous_void_fraction[i]);
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    vf_system_trans_vectors(this->void_fraction_dof_handler);
  vf_system_trans_vectors.prepare_for_coarsening_and_refinement(
    vf_set_transfer);

  // Prepare particle handle for serialization
  this->particle_handler.prepare_for_coarsening_and_refinement();

  this->pcout << "-->Repartitionning triangulation" << std::endl;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  parallel_triangulation->repartition();

  // Update cell neighbors
  if (has_periodic_boundaries)
    {
      periodic_boundaries_object.map_periodic_cells(
        *parallel_triangulation, periodic_boundaries_cells_information);

      periodic_offset =
        periodic_boundaries_object.get_periodic_offset_distance();
    }

  contact_manager.update_cell_neighbors(*parallel_triangulation,
                                        periodic_boundaries_cells_information,
                                        has_periodic_boundaries);


  boundary_cell_object.build(
    *parallel_triangulation,
    dem_parameters.floating_walls,
    dem_parameters.boundary_conditions.outlet_boundaries,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .check_for_diamond_cells,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .expand_particle_wall_contact_search,
    this->pcout);

  const auto average_minimum_maximum_cells =
    Utilities::MPI::min_max_avg(parallel_triangulation->n_active_cells(),
                                this->mpi_communicator);

  const auto average_minimum_maximum_particles = Utilities::MPI::min_max_avg(
    this->particle_handler.n_locally_owned_particles(), this->mpi_communicator);

  this->pcout << "Load balance finished " << std::endl;
  this->pcout
    << "Average, minimum and maximum number of particles on the processors are "
    << average_minimum_maximum_particles.avg << " , "
    << average_minimum_maximum_particles.min << " and "
    << average_minimum_maximum_particles.max << std::endl;
  this->pcout
    << "Minimum and maximum number of cells owned by the processors are "
    << average_minimum_maximum_cells.min << " and "
    << average_minimum_maximum_cells.max << std::endl;

  this->pcout << "Setup DOFs" << std::endl;
  this->setup_dofs();

  // Remap periodic nodes after setup of dofs
  if (has_periodic_boundaries && has_sparse_contacts)
    {
      sparse_contacts_object.map_periodic_nodes(
        this->void_fraction_constraints);
    }

  // Velocity Vectors
  std::vector<GlobalVectorType *> x_system(1 + this->previous_solutions.size());

  GlobalVectorType distributed_system(this->locally_owned_dofs,
                                      this->mpi_communicator);

  x_system[0] = &(distributed_system);

  std::vector<GlobalVectorType> distributed_previous_solutions;

  distributed_previous_solutions.reserve(this->previous_solutions.size());

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(this->locally_owned_dofs, this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }

  system_trans_vectors.interpolate(x_system);

  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      this->previous_solutions[i] = distributed_previous_solutions[i];
    }

  x_system.clear();

  // Void Fraction Vectors
  std::vector<GlobalVectorType *> vf_system(
    1 + this->previous_void_fraction.size());

  GlobalVectorType vf_distributed_system(this->locally_owned_dofs_voidfraction,
                                         this->mpi_communicator);

  vf_system[0] = &(vf_distributed_system);

  std::vector<GlobalVectorType> vf_distributed_previous_solutions;

  vf_distributed_previous_solutions.reserve(
    this->previous_void_fraction.size());

  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_distributed_previous_solutions.emplace_back(
        GlobalVectorType(this->locally_owned_dofs_voidfraction,
                         this->mpi_communicator));
      vf_system[i + 1] = &vf_distributed_previous_solutions[i];
    }

  vf_system_trans_vectors.interpolate(vf_system);

  this->nodal_void_fraction_relevant = vf_distributed_system;
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      this->previous_void_fraction[i] = vf_distributed_previous_solutions[i];
    }

  vf_system.clear();

  // Unpack particle handler after load balancing step
  this->particle_handler.unpack_after_coarsening_and_refinement();

  // Regenerate vertex to cell map
  this->vertices_cell_mapping();
}

template <int dim>
void
CFDDEMSolver<dim>::initialize_dem_parameters()
{
  this->pcout << "Initializing DEM parameters " << std::endl;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  if (has_periodic_boundaries)
    {
      periodic_boundaries_object.set_periodic_boundaries_information(
        dem_parameters.boundary_conditions.periodic_boundary_0,
        dem_parameters.boundary_conditions.periodic_boundary_1,
        dem_parameters.boundary_conditions.periodic_direction);

      periodic_boundaries_object.map_periodic_cells(
        *parallel_triangulation, periodic_boundaries_cells_information);

      // Temporary offset calculation : works only for one set of periodic
      // boundary on an axis.
      periodic_offset =
        periodic_boundaries_object.get_periodic_offset_distance();

      // Initialize the flag for particle displacement in PBC
      particle_displaced_in_pbc = false;
    }

  if (dem_parameters.model_parameters.sparse_particle_contacts)
    {
      has_sparse_contacts = true;
      sparse_contacts_object.set_parameters(
        dem_parameters.model_parameters.granular_temperature_threshold,
        dem_parameters.model_parameters.solid_fraction_threshold,
        dem_parameters.model_parameters.advect_particles);
    }

  // Finding cell neighbors
  contact_manager.execute_cell_neighbors_search(
    *parallel_triangulation,
    periodic_boundaries_cells_information,
    has_periodic_boundaries);

  // Finding boundary cells with faces
  boundary_cell_object.build(
    *parallel_triangulation,
    dem_parameters.floating_walls,
    dem_parameters.boundary_conditions.outlet_boundaries,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .check_for_diamond_cells,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .expand_particle_wall_contact_search,
    this->pcout);

  // Setting chosen contact force, insertion and integration methods
  integrator_object = set_integrator_type();
  particle_particle_contact_force_object =
    set_particle_particle_contact_force_model(
      this->cfd_dem_simulation_parameters.dem_parameters);
  particle_wall_contact_force_object = set_particle_wall_contact_force_model(
    this->cfd_dem_simulation_parameters.dem_parameters,
    *parallel_triangulation);

  this->particle_handler.sort_particles_into_subdomains_and_cells();

  displacement.resize(this->particle_handler.get_max_local_particle_index());

  force.resize(displacement.size());
  torque.resize(displacement.size());

  this->particle_handler.exchange_ghost_particles(true);

  // Updating moment of inertia container
  update_moment_of_inertia(this->particle_handler, MOI);

  this->pcout << "Finished initializing DEM parameters " << std::endl
              << "DEM time-step is " << dem_time_step << " s " << std::endl;
}

template <int dim>
void
CFDDEMSolver<dim>::update_moment_of_inertia(
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
std::shared_ptr<Integrator<dim>>
CFDDEMSolver<dim>::set_integrator_type()
{
  if (dem_parameters.model_parameters.integration_method ==
      Parameters::Lagrangian::ModelParameters::IntegrationMethod::
        velocity_verlet)
    {
      integrator_object = std::make_shared<VelocityVerletIntegrator<dim>>();
    }
  else if (dem_parameters.model_parameters.integration_method ==
           Parameters::Lagrangian::ModelParameters::IntegrationMethod::
             explicit_euler)
    {
      integrator_object = std::make_shared<ExplicitEulerIntegrator<dim>>();
    }
  else if (dem_parameters.model_parameters.integration_method ==
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
CFDDEMSolver<dim>::add_fluid_particle_interaction_force()
{
  for (auto particle = this->particle_handler.begin();
       particle != this->particle_handler.end();
       ++particle)
    {
      auto particle_properties = particle->get_properties();

      types::particle_index particle_id = particle->get_local_index();

      force[particle_id][0] +=
        particle_properties[DEM::PropertiesIndex::fem_force_x];
      force[particle_id][1] +=
        particle_properties[DEM::PropertiesIndex::fem_force_y];
      force[particle_id][2] +=
        particle_properties[DEM::PropertiesIndex::fem_force_z];
    }
}

template <int dim>
void
CFDDEMSolver<dim>::add_fluid_particle_interaction_torque()
{
  for (auto particle = this->particle_handler.begin();
       particle != this->particle_handler.end();
       ++particle)
    {
      auto particle_properties = particle->get_properties();

      types::particle_index particle_id = particle->get_local_index();

      torque[particle_id][0] +=
        particle_properties[DEM::PropertiesIndex::fem_torque_x];
      torque[particle_id][1] +=
        particle_properties[DEM::PropertiesIndex::fem_torque_y];
      torque[particle_id][2] +=
        particle_properties[DEM::PropertiesIndex::fem_torque_z];
    }
}

template <int dim>
void
CFDDEMSolver<dim>::dem_iterator(unsigned int counter)
{
  // dem_contact_build carries out the particle-particle and particle-wall
  // broad and fine searches, sort_particles_into_subdomains_and_cells, and
  // exchange_ghost
  dem_contact_build(counter);

  // Particle-particle contact force
  particle_particle_contact_force_object
    ->calculate_particle_particle_contact_force(
      contact_manager, dem_time_step, torque, force, periodic_offset);

  // Particles-walls contact force:
  particle_wall_contact_force();

  // Add fluid-particle interaction force to the force container
  add_fluid_particle_interaction_force();

  // Add fluid-particle interaction torque to the torque container
  if (this->cfd_dem_simulation_parameters.cfd_dem.rotational_viscous_torque ||
      this->cfd_dem_simulation_parameters.cfd_dem.vortical_viscous_torque)
    add_fluid_particle_interaction_torque();

  // Integration correction step (after force calculation)
  // In the first step, we have to obtain location of particles at half-step
  // time
  if (this->simulation_control->get_step_number() == 0)
    {
      integrator_object->integrate_half_step_location(
        this->particle_handler, g, dem_time_step, torque, force, MOI);
    }
  else
    {
      if (!has_sparse_contacts)
        {
          integrator_object->integrate(
            this->particle_handler, g, dem_time_step, torque, force, MOI);
        }
      else
        {
          if (counter > 0)
            {
              const auto parallel_triangulation =
                dynamic_cast<parallel::distributed::Triangulation<dim> *>(
                  &*this->triangulation);
              integrator_object->integrate(this->particle_handler,
                                           g,
                                           dem_time_step,
                                           torque,
                                           force,
                                           MOI,
                                           *parallel_triangulation,
                                           sparse_contacts_object);
            }
          else // counter == 0
            {
              // The cell average velocities and accelerations are updated
              // from the fully computed forces at step 0, but we do not use the
              // mobility status to disabled contacts at this step. The reason
              // is that the current mobility status are computed for the last
              // DEM time step (from the last CFD time step) are the force of
              // the fluid may have significantly changed the particulate
              // agitation.

              // Update the cell average velocities and accelerations
              sparse_contacts_object.update_average_velocities_acceleration(
                this->particle_handler, g, force, dem_time_step);

              integrator_object->integrate(
                this->particle_handler, g, dem_time_step, torque, force, MOI);
            }
        }
    }

  // If simulation has periodic boundaries, the particles are sorted into
  // subdomains and cells at the last DEM coupled time step otherwise the
  // particles will not match the cells that they are in when void fraction is
  // calculated with the qcm method
  if (counter == (coupling_frequency - 1))
    {
      if (has_periodic_boundaries &&
          this->cfd_dem_simulation_parameters.void_fraction->mode ==
            Parameters::VoidFractionMode::qcm)
        {
          bool particle_has_been_moved =
            periodic_boundaries_object.execute_particles_displacement(
              this->particle_handler, periodic_boundaries_cells_information);

          // Exchange information between processors
          particle_displaced_in_pbc =
            Utilities::MPI::logical_or(particle_has_been_moved,
                                       this->mpi_communicator);

          if (particle_displaced_in_pbc)
            {
              this->particle_handler.sort_particles_into_subdomains_and_cells();
              this->particle_handler.exchange_ghost_particles(true);
            }
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::dem_contact_build(unsigned int counter)
{
  // Check to see if it is contact detection step
  contact_detection_step =
    check_contact_detection_method(counter,
                                   this->cfd_dem_simulation_parameters,
                                   displacement,
                                   this->particle_handler,
                                   this->mpi_communicator,
                                   this->simulation_control,
                                   contact_detection_step,
                                   checkpoint_step,
                                   load_balance_step,
                                   smallest_contact_search_criterion);

  // Sort particles in cells
  if (contact_search_step(counter))
    {
      this->pcout << "DEM contact search at dem step " << counter << std::endl;
      contact_search_counter++;

      periodic_boundaries_object.execute_particles_displacement(
        this->particle_handler, periodic_boundaries_cells_information);

      this->particle_handler.sort_particles_into_subdomains_and_cells();

      displacement.resize(
        this->particle_handler.get_max_local_particle_index());
      force.resize(displacement.size());
      torque.resize(displacement.size());

      this->particle_handler.exchange_ghost_particles(true);

      if (has_sparse_contacts)
        {
          if (load_balance_step || checkpoint_step ||
              (this->simulation_control->is_at_start() && (counter == 0)))
            sparse_contacts_object.update_local_and_ghost_cell_set(
              this->void_fraction_dof_handler);

          sparse_contacts_object.identify_mobility_status(
            this->void_fraction_dof_handler,
            this->particle_handler,
            (*this->triangulation).n_active_cells(),
            this->mpi_communicator,
            counter);
        }

      // Updating moment of inertia container
      update_moment_of_inertia(this->particle_handler, MOI);

      // Execute broad search by filling containers of particle-particle
      // contact pair candidates and containers of particle-wall
      // contact pair candidates
      if (!(has_sparse_contacts && counter > 0))
        {
          contact_manager.execute_particle_particle_broad_search(
            this->particle_handler, has_periodic_boundaries);

          contact_manager.execute_particle_wall_broad_search(
            this->particle_handler,
            boundary_cell_object,
            solid_surfaces_mesh_info,
            dem_parameters.floating_walls,
            this->simulation_control->get_current_time());
        }
      else
        {
          contact_manager.execute_particle_particle_broad_search(
            this->particle_handler,
            sparse_contacts_object,
            has_periodic_boundaries);

          contact_manager.execute_particle_wall_broad_search(
            this->particle_handler,
            boundary_cell_object,
            solid_surfaces_mesh_info,
            dem_parameters.floating_walls,
            this->simulation_control->get_current_time(),
            sparse_contacts_object);
        }

      // Update contacts, remove replicates and add new contact pairs
      // to the contact containers when particles are exchanged between
      // processors
      contact_manager.update_contacts(has_periodic_boundaries);

      // Updates the iterators to particles in local-local contact
      // containers
      contact_manager.update_local_particles_in_cells(this->particle_handler,
                                                      load_balance_step,
                                                      has_periodic_boundaries);

      // Execute fine search by updating particle-particle contact
      // containers regards the neighborhood threshold
      contact_manager.execute_particle_particle_fine_search(
        neighborhood_threshold_squared,
        has_periodic_boundaries,
        periodic_offset);

      // Execute fine search by updating particle-wall contact containers
      // regards the neighborhood threshold
      contact_manager.execute_particle_wall_fine_search(
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        neighborhood_threshold_squared);

      // Reset different steps. The contact build should be performed everytime
      // we restart the simulation or everytime load balancing is performed. At
      // the end of the restart step or the load balance step, and after all
      // necessary contact build, vector resizing, and solution transfers have
      // been performed, this functions are set to false. The checkpoint_step
      // remains false for the duration of the simulation while the load
      // balancing is reset everytime load balancing is called.
      checkpoint_step           = false;
      load_balance_step         = false;
      particle_displaced_in_pbc = false;
    }
  else
    {
      this->particle_handler.update_ghost_particles();
    }
  // TODO add DEM post-processing
}

template <int dim>
void
CFDDEMSolver<dim>::write_DEM_output_results()
{
  const std::string folder = dem_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    dem_parameters.simulation_control.output_name + "_particles";
  const unsigned int iter = this->simulation_control->get_step_number();
  const double       time = this->simulation_control->get_current_time();
  const unsigned int group_files =
    dem_parameters.simulation_control.group_files;

  // Write particles
  Visualization<dim> particle_data_out;
  particle_data_out.build_patches(this->particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler,
                            particle_data_out,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            this->mpi_communicator);
}

template <int dim>
void
CFDDEMSolver<dim>::particle_wall_contact_force()
{
  // Particle-wall contact force
  particle_wall_contact_force_object->calculate_particle_wall_contact_force(
    contact_manager.particle_wall_in_contact, dem_time_step, torque, force);

  if (this->cfd_dem_simulation_parameters.dem_parameters.forces_torques
        .calculate_force_torque)
    {
      forces_boundary_information[this->simulation_control->get_step_number()] =
        particle_wall_contact_force_object->get_force();
      torques_boundary_information[this->simulation_control
                                     ->get_step_number()] =
        particle_wall_contact_force_object->get_torque();
    }

  // Particle-floating wall contact force
  if (dem_parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_contact_force_object->calculate_particle_wall_contact_force(
        contact_manager.particle_floating_wall_in_contact,
        dem_time_step,
        torque,
        force);
    }

  particle_point_line_contact_force_object
    .calculate_particle_point_contact_force(
      &contact_manager.particle_points_in_contact,
      dem_parameters.lagrangian_physical_properties,
      force);

  if constexpr (dim == 3)
    {
      particle_point_line_contact_force_object
        .calculate_particle_line_contact_force(
          &contact_manager.particle_lines_in_contact,
          dem_parameters.lagrangian_physical_properties,
          force);
    }
}


template <int dim>
void
CFDDEMSolver<dim>::dem_post_process_results()
{
  // Update statistics on contact list
  double number_of_list_built_since_last_log =
    double(contact_search_total_number) - contact_list.total;
  contact_list.max =
    std::max(number_of_list_built_since_last_log, contact_list.max);
  contact_list.min =
    std::min(number_of_list_built_since_last_log, contact_list.min);
  contact_list.total   = double(contact_search_total_number);
  contact_list.average = contact_list.total /
                         (this->simulation_control->get_step_number()) *
                         this->simulation_control->get_log_frequency();

  // Calculate statistics on the particles
  statistics translational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::dem_statistic_variable::translational_kinetic_energy>(
    this->particle_handler, this->mpi_communicator);
  statistics rotational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::dem_statistic_variable::rotational_kinetic_energy>(
    this->particle_handler, this->mpi_communicator);
  statistics velocity =
    calculate_granular_statistics<dim, DEM::dem_statistic_variable::velocity>(
      this->particle_handler, this->mpi_communicator);
  statistics omega =
    calculate_granular_statistics<dim, DEM::dem_statistic_variable::omega>(
      this->particle_handler, this->mpi_communicator);

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

      announce_string(this->pcout, "Particle statistics");
      report.write_text(std::cout, dealii::TableHandler::org_mode_table);
    }

  if (dem_parameters.post_processing.Lagrangian_post_processing &&
      this->simulation_control->is_output_iteration())
    {
      const auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation);


      dem_post_processing_object.write_post_processing_results(
        *parallel_triangulation,
        grid_pvdhandler,
        this->dof_handler,
        this->particle_handler,
        dem_parameters,
        this->simulation_control->get_current_time(),
        this->simulation_control->get_step_number(),
        this->mpi_communicator,
        sparse_contacts_object);
    }
}

template <int dim>
void
CFDDEMSolver<dim>::postprocess_fd(bool first_iteration)
{
  this->pcout
    << "---------------------------------------------------------------"
    << std::endl;

  this->GLSNavierStokesSolver<dim>::postprocess_fd(first_iteration);

  // Visualization
  if (this->simulation_control->is_output_iteration())
    {
      write_DEM_output_results();
    }
}

template <int dim>
void
CFDDEMSolver<dim>::postprocess_cfd_dem()
{
  // Calculate total volume of fluid and solid
  if (this->simulation_parameters.post_processing.calculate_phase_volumes)
    {
      TimerOutput::Scope t(this->computing_timer, "total_volume_calculation");
      double             total_volume_fluid, total_volume_particles;
      std::tie(total_volume_fluid, total_volume_particles) =
        calculate_fluid_and_particle_volumes(this->void_fraction_dof_handler,
                                             this->nodal_void_fraction_relevant,
                                             *this->cell_quadrature,
                                             *this->mapping);
      this->table_phase_volumes.add_value(
        "time", this->simulation_control->get_current_time());
      this->table_phase_volumes.add_value("total-volume-fluid",
                                          total_volume_fluid);
      this->table_phase_volumes.add_value("total-volume-particles",
                                          total_volume_particles);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Total volume of fluid: "
                      << std::setprecision(
                           this->simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           total_volume_fluid
                      << " m^3" << std::endl;
          this->pcout << "Total volume of particles: "
                      << std::setprecision(
                           this->simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           total_volume_particles
                      << " m^3" << std::endl;
        }

      // Output pressure drop to a text file from processor 0
      if ((this->simulation_control->get_step_number() %
             this->simulation_parameters.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            this->simulation_parameters.simulation_control.output_folder +
            this->simulation_parameters.post_processing
              .phase_volumes_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          table_phase_volumes.set_precision("time", 12);
          table_phase_volumes.set_precision("total-volume-fluid", 12);
          table_phase_volumes.set_precision("total-volume-particles", 12);
          this->table_phase_volumes.write_text(output);
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::dynamic_flow_control()
{
  if (this->simulation_parameters.flow_control.enable_flow_control &&
      this->simulation_parameters.simulation_control.method !=
        Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      // Calculate the average velocity according to the void fraction
      unsigned int flow_direction =
        this->simulation_parameters.flow_control.flow_direction;
      double average_velocity =
        calculate_average_velocity(this->dof_handler,
                                   this->void_fraction_dof_handler,
                                   this->present_solution,
                                   this->nodal_void_fraction_relevant,
                                   flow_direction,
                                   *this->cell_quadrature,
                                   *this->mapping);

      // Calculate the beta force for fluid
      this->flow_control.calculate_beta(
        average_velocity,
        this->simulation_control->get_time_step(),
        this->simulation_control->get_step_number());

      Tensor<1, 3> beta_particle;
      if (this->simulation_parameters.flow_control.enable_beta_particle)
        {
          // Calculate the beta for particles and add it to force tensor g
          AssertThrow(
            dem_parameters.lagrangian_physical_properties
                .particle_type_number == 1,
            ExcMessage(std::string(
              "Beta force calculation for particles is not implemented for multiple particle types.")));

          double fluid_density =
            this->simulation_parameters.physical_properties_manager
              .get_density_scale();
          double particle_density =
            dem_parameters.lagrangian_physical_properties.density_particle[0];
          beta_particle =
            this->flow_control.get_beta_particles(fluid_density,
                                                  particle_density);
          g = dem_parameters.lagrangian_physical_properties.g + beta_particle;
        }

      // Showing results
      if (this->simulation_parameters.flow_control.verbosity ==
            Parameters::Verbosity::verbose &&
          this->simulation_control->get_step_number() > 0 &&
          this->this_mpi_process == 0)
        {
          announce_string(this->pcout, "Flow control summary");
          this->pcout << "Fluid space-average velocity: " << average_velocity
                      << std::endl;
          this->pcout << "Fluid beta force: "
                      << this->flow_control.get_beta()[flow_direction]
                      << std::endl;
          this->pcout << "Particle beta force: "
                      << beta_particle[flow_direction] << std::endl;
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::print_particles_summary()
{
  this->pcout << "Particle Summary" << std::endl;
  this->pcout << "id, x, y, z, v_x, v_y, v_z " << std::endl;
  std::map<int, Particles::ParticleIterator<dim>> global_particles;
  unsigned int                                    current_id, id_max = 0;

  // Mapping of all particles & find the max id on current processor
  for (auto particle = this->particle_handler.begin();
       particle != this->particle_handler.end();
       ++particle)
    {
      current_id = particle->get_id();
      id_max     = std::max(current_id, id_max);

      global_particles.insert({current_id, particle});
    }

  for (unsigned int i = 0; i <= id_max; i++)
    {
      for (auto &iterator : global_particles)
        {
          unsigned int id = iterator.first;
          if (id == i)
            {
              auto particle            = iterator.second;
              auto particle_properties = particle->get_properties();
              auto particle_location   = particle->get_location();

              std::stringstream ss;
              ss << std::setprecision(6) << id << " " << particle_location
                 << " " << particle_properties[DEM::PropertiesIndex::v_x] << " "
                 << particle_properties[DEM::PropertiesIndex::v_y] << " "
                 << particle_properties[DEM::PropertiesIndex::v_z];
              this->pcout << ss.str() << std::endl;
            }
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::dem_setup_contact_parameters()
{
  coupling_frequency =
    this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency;

  // Initialize DEM Parameters
  dem_parameters.lagrangian_physical_properties =
    this->cfd_dem_simulation_parameters.dem_parameters
      .lagrangian_physical_properties;
  g = dem_parameters.lagrangian_physical_properties.g;
  dem_parameters.boundary_conditions =
    this->cfd_dem_simulation_parameters.dem_parameters.boundary_conditions;
  dem_parameters.floating_walls =
    this->cfd_dem_simulation_parameters.dem_parameters.floating_walls;
  dem_parameters.model_parameters =
    this->cfd_dem_simulation_parameters.dem_parameters.model_parameters;
  dem_parameters.simulation_control =
    this->cfd_dem_simulation_parameters.dem_parameters.simulation_control;
  dem_parameters.post_processing =
    this->cfd_dem_simulation_parameters.dem_parameters.post_processing;
  dem_parameters.mesh = this->cfd_dem_simulation_parameters.dem_parameters.mesh;
  dem_parameters.restart =
    this->cfd_dem_simulation_parameters.dem_parameters.restart;
  size_distribution_object_container.reserve(
    dem_parameters.lagrangian_physical_properties.particle_type_number);

  maximum_particle_diameter = 0.;
  for (unsigned int particle_type = 0;
       particle_type <
       dem_parameters.lagrangian_physical_properties.particle_type_number;
       particle_type++)
    {
      if (dem_parameters.lagrangian_physical_properties.distribution_type.at(
            particle_type) ==
          Parameters::Lagrangian::SizeDistributionType::uniform)
        {
          size_distribution_object_container.push_back(
            std::make_shared<UniformDistribution>(
              dem_parameters.lagrangian_physical_properties
                .particle_average_diameter.at(particle_type)));
        }
      else if (dem_parameters.lagrangian_physical_properties.distribution_type
                 .at(particle_type) ==
               Parameters::Lagrangian::SizeDistributionType::normal)
        {
          size_distribution_object_container.push_back(
            std::make_shared<NormalDistribution>(
              dem_parameters.lagrangian_physical_properties
                .particle_average_diameter.at(particle_type),
              dem_parameters.lagrangian_physical_properties.particle_size_std
                .at(particle_type),
              dem_parameters.lagrangian_physical_properties
                .seed_for_distributions[particle_type]));
        }
      else if (dem_parameters.lagrangian_physical_properties.distribution_type
                 .at(particle_type) ==
               Parameters::Lagrangian::SizeDistributionType::custom)
        {
          size_distribution_object_container.push_back(
            std::make_shared<CustomDistribution>(
              dem_parameters.lagrangian_physical_properties
                .particle_custom_diameter.at(particle_type),
              dem_parameters.lagrangian_physical_properties
                .particle_custom_probability.at(particle_type),
              dem_parameters.lagrangian_physical_properties
                .seed_for_distributions[particle_type]));
        }

      maximum_particle_diameter = std::max(
        maximum_particle_diameter,
        size_distribution_object_container[particle_type]->find_max_diameter());
    }

  neighborhood_threshold_squared =
    std::pow(dem_parameters.model_parameters.neighborhood_threshold *
               maximum_particle_diameter,
             2);


  // Finding the smallest contact search frequency criterion between (smallest
  // cell size - largest particle radius) and (security factor * (blob diameter
  // - 1) *  largest particle radius). This value is used in
  // find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(*this->triangulation) -
              maximum_particle_diameter * 0.5),
             (dem_parameters.model_parameters.dynamic_contact_search_factor *
              (dem_parameters.model_parameters.neighborhood_threshold - 1) *
              maximum_particle_diameter * 0.5));

  dem_time_step =
    this->simulation_control->get_time_step() / coupling_frequency;

  double rayleigh_time_step = 0;

  for (unsigned int i = 0;
       i < dem_parameters.lagrangian_physical_properties.particle_type_number;
       ++i)
    rayleigh_time_step = std::max(
      M_PI_2 *
        dem_parameters.lagrangian_physical_properties
          .particle_average_diameter[i] *
        sqrt(2 *
             dem_parameters.lagrangian_physical_properties.density_particle[i] *
             (2 + dem_parameters.lagrangian_physical_properties
                    .poisson_ratio_particle[i]) *
             (1 - dem_parameters.lagrangian_physical_properties
                    .poisson_ratio_particle[i]) /
             dem_parameters.lagrangian_physical_properties
               .youngs_modulus_particle[i]) /
        (0.1631 * dem_parameters.lagrangian_physical_properties
                    .poisson_ratio_particle[i] +
         0.8766),
      rayleigh_time_step);

  const double time_step_rayleigh_ratio = dem_time_step / rayleigh_time_step;
  this->pcout << "DEM time-step is " << time_step_rayleigh_ratio * 100
              << "% of Rayleigh time step" << std::endl;

  // Initialize contact detection step
  contact_detection_step = false;
  load_balance_step      = false;

  // Check if there are periodic boundaries
  for (unsigned int i_bc = 0;
       i_bc < dem_parameters.boundary_conditions.bc_types.size();
       ++i_bc)
    {
      if (dem_parameters.boundary_conditions.bc_types[i_bc] ==
          Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
        {
          has_periodic_boundaries = true;
          break;
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::manage_triangulation_connections()
{
  // Necessary signals for load balancing. This signals are only connected if
  // load balancing is enabled. This helps prevent errors in read_dem function
  // for Deal.II version 9.3
  //  const auto parallel_triangulation =
  //    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
  //      &*this->triangulation);

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  if (dem_parameters.model_parameters.load_balance_method !=
      Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::none)
    {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 6)
      parallel_triangulation->signals.weight.connect(
        [](const typename Triangulation<dim>::cell_iterator &,
           const typename Triangulation<dim>::CellStatus) -> unsigned int {
          return 1000;
        });

      parallel_triangulation->signals.weight.connect(
        [&](const typename parallel::distributed::Triangulation<
              dim>::cell_iterator &cell,
            const typename parallel::distributed::Triangulation<dim>::CellStatus
              status) -> unsigned int {
          return this->cell_weight(cell, status);
        });

#else
      parallel_triangulation->signals.weight.connect(
        [](const typename Triangulation<dim>::cell_iterator &,
           const CellStatus) -> unsigned int { return 1000; });

      parallel_triangulation->signals.weight.connect(
        [&](const typename parallel::distributed::Triangulation<
              dim>::cell_iterator &cell,
            const CellStatus       status) -> unsigned int {
          return this->cell_weight(cell, status);
        });
#endif
    }
}

template <int dim>
void
CFDDEMSolver<dim>::solve()
{
  read_mesh_and_manifolds(
    *this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  manage_triangulation_connections();

  dem_setup_contact_parameters();

  // Reading DEM start file information
  if (this->cfd_dem_simulation_parameters.void_fraction->read_dem == true &&
      this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
          .restart == false)
    read_dem();

  this->setup_dofs();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  // In the case the simulation is being restarted from a checkpoint file, the
  // checkpoint_step parameter is set to true. This allows to perform all
  // operations related to restarting a simulation. Once all operations have
  // been performed, this checkpoint_step is reset to false. It is only set once
  // and reset once since restarting only occurs once.
  if (this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart == true)
    checkpoint_step = true;
  else
    checkpoint_step = false;

  // Initialize the DEM parameters and generate the required ghost particles
  initialize_dem_parameters();
  this->particle_handler.exchange_ghost_particles(true);


  // Remap periodic nodes after setup of dofs
  if (has_periodic_boundaries && has_sparse_contacts)
    {
      sparse_contacts_object.map_periodic_nodes(
        this->void_fraction_constraints);
    }

  // Calculate first instance of void fraction once particles are set-up
  this->vertices_cell_mapping();
  if (!checkpoint_step)
    this->initialize_void_fraction();

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      bool refinement_step;
      if (this->simulation_parameters.mesh_adaptation.refinement_at_frequency)
        refinement_step =
          this->simulation_control->get_step_number() %
            this->simulation_parameters.mesh_adaptation.frequency !=
          0;
      else
        refinement_step = this->simulation_control->get_step_number() == 0;
      if (refinement_step ||
          this->simulation_parameters.mesh_adaptation.type ==
            Parameters::MeshAdaptation::Type::none ||
          this->simulation_control->is_at_start())
        {
          // We allow the physics to update their boundary conditions
          // according to their own parameters
          this->update_boundary_conditions();
          this->multiphysics->update_boundary_conditions();
        }

      this->dynamic_flow_control();

      if (!this->simulation_control->is_at_start())
        {
          NavierStokesBase<dim, GlobalVectorType, IndexSet>::refine_mesh();
          this->vertices_cell_mapping();
        }

      this->calculate_void_fraction(
        this->simulation_control->get_current_time(), load_balance_step);
      this->iterate();

      if (this->cfd_dem_simulation_parameters.cfd_parameters.test.enabled)
        {
          switch (
            this->cfd_dem_simulation_parameters.cfd_parameters.test.test_type)
            {
              case Parameters::Testing::TestType::particles:
                {
                  print_particles_summary();
                  break;
                }
              case Parameters::Testing::TestType::mobility_status:
                {
                  if (this->simulation_control->is_at_end())

                    {
                      // Get mobility status vector sorted by cell id
                      Vector<float> mobility_status(
                        this->triangulation->n_active_cells());
                      sparse_contacts_object.get_mobility_status_vector(
                        mobility_status);

                      // Output mobility status vector
                      visualization_object.print_intermediate_format(
                        mobility_status,
                        this->void_fraction_dof_handler,
                        this->mpi_communicator);
                    }
                  break;
                }
              default:
                print_particles_summary();
            }
        }

      {
        announce_string(this->pcout, "DEM");
        TimerOutput::Scope t(this->computing_timer, "DEM_Iterator");

        contact_search_counter = 0;
        for (unsigned int dem_counter = 0; dem_counter < coupling_frequency;
             ++dem_counter)
          {
            // dem_iterator carries out the particle-particle and
            // particle_wall force calculations, integration and
            // update_ghost
            dem_iterator(dem_counter);
          }
      }
      contact_search_total_number += contact_search_counter;

      this->pcout << "Finished " << coupling_frequency << " DEM iterations "
                  << std::endl;

      this->postprocess(false);
      this->postprocess_cfd_dem();
      this->finish_time_step_fd();

      this->GLSVANSSolver<dim>::monitor_mass_conservation();

      if (this->cfd_dem_simulation_parameters.cfd_dem.particle_statistics)
        dem_post_process_results();

      // Load balancing
      // The input argument to this function is set to zero as this integer is
      // not used for the check_load_balance_step function and is only important
      // for the check_contact_search_step function.
      load_balance_step =
        check_load_balance_method(this->cfd_dem_simulation_parameters,
                                  this->particle_handler,
                                  this->mpi_communicator,
                                  this->n_mpi_processes,
                                  this->simulation_control);

      if (load_balance_step || checkpoint_step)
        {
          load_balance();
        }
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D CFD-DEM solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class CFDDEMSolver<2>;
template class CFDDEMSolver<3>;
