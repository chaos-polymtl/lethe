#include <core/solutions_output.h>

#include <dem/dem_solver_parameters.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/find_maximum_particle_size.h>
#include <dem/gear3_integrator.h>
#include <dem/pp_linear_force.h>
#include <dem/pp_nonlinear_force.h>
#include <dem/pw_linear_force.h>
#include <dem/pw_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>
#include <fem-dem/cfd_dem_coupling.h>

// Constructor for class CFD-DEM
template <int dim>
CFDDEMSolver<dim>::CFDDEMSolver(CFDDEMSimulationParameters<dim> &nsparam)
  : GLSVANSSolver<dim>(nsparam)
{
  coupling_frequency =
    this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency;

  contact_detection_frequency =
    this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
      .contact_detection_frequency;

  standard_deviation_multiplier = 2.5;

  load_balancing_frequency = this->cfd_dem_simulation_parameters.dem_parameters
                               .model_parameters.load_balance_frequency;

  if (this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart == true)
    checkpoint_step = true;

  maximum_particle_diameter =
    find_maximum_particle_size(this->cfd_dem_simulation_parameters
                                 .dem_parameters.lagrangian_physical_properties,
                               standard_deviation_multiplier);
  neighborhood_threshold_squared =
    std::pow(this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
                 .neighborhood_threshold *
               maximum_particle_diameter,
             2);

  if (this->cfd_dem_simulation_parameters.dem_parameters.mesh.type ==
      Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        tria,
        this->cfd_dem_simulation_parameters.dem_parameters.mesh.grid_type,
        this->cfd_dem_simulation_parameters.dem_parameters.mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  triangulation_cell_diameter = 0.5 * GridTools::diameter(tria);

  // Finding the smallest contact search frequency criterion between (smallest
  // cell size - largest particle radius) and (security factor * (blab
  // diamater
  // - 1) *  largest particle radius). This value is used in
  // find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(tria) -
              maximum_particle_diameter * 0.5),
             (this->cfd_dem_simulation_parameters.dem_parameters
                .model_parameters.dynamic_contact_search_factor *
              (this->cfd_dem_simulation_parameters.dem_parameters
                 .model_parameters.neighborhood_threshold -
               1) *
              maximum_particle_diameter * 0.5));

  dem_parameters.lagrangian_physical_properties =
    this->cfd_dem_simulation_parameters.dem_parameters
      .lagrangian_physical_properties;
  dem_parameters.boundary_conditions =
    this->cfd_dem_simulation_parameters.dem_parameters.boundary_conditions;
  dem_parameters.floating_walls =
    this->cfd_dem_simulation_parameters.dem_parameters.floating_walls;
  dem_parameters.model_parameters =
    this->cfd_dem_simulation_parameters.dem_parameters.model_parameters;

  dem_time_step =
    this->simulation_control->get_time_step() / coupling_frequency;

  if (dem_parameters.model_parameters.contact_detection_method ==
      Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::constant)
    {
      check_contact_search_step =
        &CFDDEMSolver<dim>::check_contact_search_step_constant;
    }
  else if (dem_parameters.model_parameters.contact_detection_method ==
           Parameters::Lagrangian::ModelParameters::ContactDetectionMethod::
             dynamic)
    {
      check_contact_search_step =
        &CFDDEMSolver<dim>::check_contact_search_step_dynamic;
    }
  else
    {
      throw std::runtime_error(
        "Specified contact detection method is not valid");
    }

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
  // These connections only need to be created once, so we
  // have set them up in the constructor of this class.
  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  parallel_triangulation->signals.cell_weight.connect(
    [&](const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int { return this->cell_weight(cell, status); });

  parallel_triangulation->signals.pre_distributed_repartition.connect(std::bind(
    &Particles::ParticleHandler<dim>::register_store_callback_function,
    &this->particle_handler));

  parallel_triangulation->signals.post_distributed_repartition.connect(
    std::bind(&Particles::ParticleHandler<dim>::register_load_callback_function,
              &this->particle_handler,
              false));
  parallel_triangulation->signals.pre_distributed_refinement.connect(std::bind(
    &Particles::ParticleHandler<dim>::register_store_callback_function,
    &this->particle_handler));

  parallel_triangulation->signals.post_distributed_refinement.connect(
    std::bind(&Particles::ParticleHandler<dim>::register_load_callback_function,
              &this->particle_handler,
              false));

  // Necessary signals for writing and reading checkpoints
  parallel_triangulation->signals.pre_distributed_save.connect(std::bind(
    &Particles::ParticleHandler<dim>::register_store_callback_function,
    &this->particle_handler));

  parallel_triangulation->signals.post_distributed_load.connect(
    std::bind(&Particles::ParticleHandler<dim>::register_load_callback_function,
              &this->particle_handler,
              true));

  // Initilize contact detection step
  contact_detection_step = true;
  checkpoint_step        = false;
  load_balance_step      = false;

  // Setting load-balance method (single-step, frequent or dynamic)
  if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
        .load_balance_method ==
      Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::once)
    {
      check_load_balance_step = &CFDDEMSolver<dim>::check_load_balance_once;
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::frequent)
    {
      check_load_balance_step = &CFDDEMSolver<dim>::check_load_balance_frequent;
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::dynamic)
    {
      check_load_balance_step = &CFDDEMSolver<dim>::check_load_balance_dynamic;
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .load_balance_method ==
           Parameters::Lagrangian::ModelParameters::LoadBalanceMethod::none)
    {
      check_load_balance_step = &CFDDEMSolver<dim>::no_load_balance;
    }
  else
    {
      throw std::runtime_error("Specified load balance method is not valid");
    }
}

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

  this->pcout << "Finished reading DEM checkpoint " << std::endl
              << this->particle_handler.n_global_particles()
              << " particles are in the simulation " << std::endl;

  write_DEM_output_results();
}

template <int dim>
void
CFDDEMSolver<dim>::write_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "Write_Checkpoint");
  this->pcout << "Writing restart file" << std::endl;

  std::string prefix = this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);
      this->pvdhandler.save(prefix);
      particles_pvdhandler.save(prefix_particles);
    }

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << this->particle_handler;

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;

  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;
  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<const TrilinosWrappers::MPI::Vector *> av_set_transfer =
        this->average_velocities->save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    system_trans_vectors(this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);


  // Void Fraction
  std::vector<const TrilinosWrappers::MPI::Vector *> vf_set_transfer;
  vf_set_transfer.push_back(&this->nodal_void_fraction_relevant);
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_set_transfer.push_back(&this->previous_void_fraction[i]);
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    vf_system_trans_vectors(this->void_fraction_dof_handler);
  vf_system_trans_vectors.prepare_for_serialization(vf_set_transfer);

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      std::string triangulationName = prefix + ".triangulation";
      parallel_triangulation->save(prefix + ".triangulation");
    }

  this->multiphysics->write_checkpoint();
}

template <int dim>
void
CFDDEMSolver<dim>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "read_checkpoint");
  std::string prefix = this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";

  this->simulation_control->read(prefix);
  this->pvdhandler.read(prefix);
  particles_pvdhandler.read(prefix_particles);

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

  // Velocity Vectors
  std::vector<TrilinosWrappers::MPI::Vector *> x_system(
    1 + this->previous_solutions.size());

  TrilinosWrappers::MPI::Vector distributed_system(this->locally_owned_dofs,
                                                   this->mpi_communicator);

  x_system[0] = &(distributed_system);

  std::vector<TrilinosWrappers::MPI::Vector> distributed_previous_solutions;

  distributed_previous_solutions.reserve(this->previous_solutions.size());

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(this->locally_owned_dofs,
                                      this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }


  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    system_trans_vectors(this->dof_handler);

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<TrilinosWrappers::MPI::Vector *> sum_vectors =
        this->average_velocities->read(prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);

  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      this->previous_solutions[i] = distributed_previous_solutions[i];
    }

  x_system.clear();

  // Void Fraction Vectors
  std::vector<TrilinosWrappers::MPI::Vector *> vf_system(
    1 + this->previous_void_fraction.size());

  TrilinosWrappers::MPI::Vector vf_distributed_system(
    this->locally_owned_dofs_voidfraction, this->mpi_communicator);

  vf_system[0] = &(vf_distributed_system);

  std::vector<TrilinosWrappers::MPI::Vector> vf_distributed_previous_solutions;

  vf_distributed_previous_solutions.reserve(
    this->previous_void_fraction.size());

  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(this->locally_owned_dofs_voidfraction,
                                      this->mpi_communicator));
      vf_system[i + 1] = &vf_distributed_previous_solutions[i];
    }

  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    vf_system_trans_vectors(this->void_fraction_dof_handler);

  vf_system_trans_vectors.deserialize(vf_system);

  this->nodal_void_fraction_relevant = vf_distributed_system;
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      this->previous_void_fraction[i] = vf_distributed_previous_solutions[i];
    }

  vf_system.clear();

  this->multiphysics->read_checkpoint();
}

template <int dim>
unsigned int
CFDDEMSolver<dim>::cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const typename parallel::distributed::Triangulation<dim>::CellStatus status)
  const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // This determines how important particle work is compared to cell
  // work (by default every cell has a weight of 1000).
  // We set the weight per particle much higher to indicate that
  // the particle load is more important than the fluid load. The optimal
  // value of this number depends on the application and can range from 0
  // (cheap particle operations, expensive cell operations) to much larger
  // than 1000 (expensive particle operations, cheap cell operations, like in
  // this case). This parameter will need to be tuned for different cases of
  // CFD-DEM coupling.
  const unsigned int particle_weight = 2000;

  // This does not use adaptive refinement, therefore every cell
  // should have the status CELL_PERSIST. However this function can also
  // be used to distribute load during refinement, therefore we consider
  // refined or coarsened cells as well.
  if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
      status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
    {
      const unsigned int n_particles_in_cell =
        this->particle_handler.n_particles_in_cell(cell);
      return n_particles_in_cell * particle_weight;
    }
  else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
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
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;
  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    system_trans_vectors(this->dof_handler);
  system_trans_vectors.prepare_for_coarsening_and_refinement(sol_set_transfer);


  // Void Fraction
  std::vector<const TrilinosWrappers::MPI::Vector *> vf_set_transfer;
  vf_set_transfer.push_back(&this->nodal_void_fraction_relevant);
  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_set_transfer.push_back(&this->previous_void_fraction[i]);
    }

  // Prepare for Serialization
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    vf_system_trans_vectors(this->void_fraction_dof_handler);
  vf_system_trans_vectors.prepare_for_coarsening_and_refinement(
    vf_set_transfer);


  this->pcout << "-->Repartitionning triangulation" << std::endl;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  parallel_triangulation->repartition();

  cells_local_neighbor_list.clear();
  cells_ghost_neighbor_list.clear();

  cell_neighbors_object.find_cell_neighbors(*parallel_triangulation,
                                            cells_local_neighbor_list,
                                            cells_ghost_neighbor_list);

  boundary_cell_object.build(
    *parallel_triangulation,
    dem_parameters.floating_walls,
    dem_parameters.boundary_conditions.outlet_boundaries,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .check_for_diamond_cells,
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

  // Velocity Vectors
  std::vector<TrilinosWrappers::MPI::Vector *> x_system(
    1 + this->previous_solutions.size());

  TrilinosWrappers::MPI::Vector distributed_system(this->locally_owned_dofs,
                                                   this->mpi_communicator);

  x_system[0] = &(distributed_system);

  std::vector<TrilinosWrappers::MPI::Vector> distributed_previous_solutions;

  distributed_previous_solutions.reserve(this->previous_solutions.size());

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(this->locally_owned_dofs,
                                      this->mpi_communicator));
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
  std::vector<TrilinosWrappers::MPI::Vector *> vf_system(
    1 + this->previous_void_fraction.size());

  TrilinosWrappers::MPI::Vector vf_distributed_system(
    this->locally_owned_dofs_voidfraction, this->mpi_communicator);

  vf_system[0] = &(vf_distributed_system);

  std::vector<TrilinosWrappers::MPI::Vector> vf_distributed_previous_solutions;

  vf_distributed_previous_solutions.reserve(
    this->previous_void_fraction.size());

  for (unsigned int i = 0; i < this->previous_void_fraction.size(); ++i)
    {
      vf_distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(this->locally_owned_dofs_voidfraction,
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
}

template <int dim>
inline bool
CFDDEMSolver<dim>::no_load_balance(const unsigned int &)
{
  return false;
}

template <int dim>
inline bool
CFDDEMSolver<dim>::check_contact_search_step_constant(
  const unsigned int &counter)
{
  return ((counter % contact_detection_frequency) == 0);
}

template <int dim>
inline bool
CFDDEMSolver<dim>::check_contact_search_step_dynamic(const unsigned int &)
{
  // The sorting into subdomain step checks whether or not the current time step
  // is a step that requires sorting particles into subdomains and cells. This
  // is applicable if any of the following three conditions apply:if its a load
  // balancing step, a restart simulation step, or a contact detection tsep.
  bool sorting_in_subdomains_step =
    (checkpoint_step || load_balance_step || contact_detection_step);

  contact_detection_step = find_contact_detection_step<dim>(
    this->particle_handler,
    this->simulation_control->get_time_step() /
      this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency,
    smallest_contact_search_criterion,
    this->mpi_communicator,
    sorting_in_subdomains_step,
    displacement);

  return contact_detection_step;
}


template <int dim>
inline bool
CFDDEMSolver<dim>::check_load_balance_once(const unsigned int &)
{
  bool load_balance_step = (this->simulation_control->get_step_number() ==
                            this->cfd_dem_simulation_parameters.dem_parameters
                              .model_parameters.load_balance_step);

  if (load_balance_step || checkpoint_step)
    {
      load_balance();
    }

  return load_balance_step;
}

template <int dim>
inline bool
CFDDEMSolver<dim>::check_load_balance_frequent(const unsigned int &)
{
  bool load_balance_step = (this->simulation_control->get_step_number() %
                              this->cfd_dem_simulation_parameters.dem_parameters
                                .model_parameters.load_balance_frequency ==
                            0);

  if (load_balance_step || checkpoint_step)
    {
      load_balance();
    }

  return load_balance_step;
}

template <int dim>
inline bool
CFDDEMSolver<dim>::check_load_balance_dynamic(const unsigned int &)
{
  bool load_balance_step = false;
  if (this->simulation_control->get_step_number() %
        this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
          .dynamic_load_balance_check_frequency ==
      0)
    {
      unsigned int maximum_particle_number_on_proc = 0;
      unsigned int minimum_particle_number_on_proc = 0;

      maximum_particle_number_on_proc =
        Utilities::MPI::max(this->particle_handler.n_locally_owned_particles(),
                            this->mpi_communicator);
      minimum_particle_number_on_proc =
        Utilities::MPI::min(this->particle_handler.n_locally_owned_particles(),
                            this->mpi_communicator);

      if ((maximum_particle_number_on_proc - minimum_particle_number_on_proc) >
            this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
                .load_balance_threshold *
              (this->particle_handler.n_global_particles() /
               this->n_mpi_processes) ||
          checkpoint_step)
        {
          load_balance();
          load_balance_step = true;
        }
    }

  return load_balance_step;
}

template <int dim>
void
CFDDEMSolver<dim>::initialize_dem_parameters()
{
  this->pcout << "Initializing DEM parameters " << std::endl;

  // TODO write read checkpoint for CFD-DEM

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  // Finding cell neighbors
  cell_neighbors_object.find_cell_neighbors(*parallel_triangulation,
                                            cells_local_neighbor_list,
                                            cells_ghost_neighbor_list);
  // Finding boundary cells with faces
  boundary_cell_object.build(
    *parallel_triangulation,
    dem_parameters.floating_walls,
    dem_parameters.boundary_conditions.outlet_boundaries,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .check_for_diamond_cells,
    this->pcout);

  // Setting chosen contact force, insertion and integration methods
  integrator_object       = set_integrator_type();
  pp_contact_force_object = set_pp_contact_force();
  pw_contact_force_object = set_pw_contact_force();

  this->particle_handler.sort_particles_into_subdomains_and_cells();

#if DEAL_II_VERSION_GTE(10, 0, 0)
  displacement.resize(this->particle_handler.get_max_local_particle_index());
#else
  {
    unsigned int max_particle_id = 0;
    for (const auto &particle : this->particle_handler)
      max_particle_id = std::max(max_particle_id, particle.get_id());
    displacement.resize(max_particle_id + 1);
  }
#endif

  force.resize(displacement.size());
  momentum.resize(displacement.size());


  this->particle_handler.exchange_ghost_particles(true);

  // Updating moment of inertia container
  update_moment_of_inertia(this->particle_handler, MOI);

  this->pcout << "Finished initializing DEM parameters " << std::endl
              << "DEM time-step is " << dem_time_step << " s ";
}

template <int dim>
void
CFDDEMSolver<dim>::update_moment_of_inertia(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  std::vector<double> &                    MOI)
{
  MOI.resize(momentum.size());

  for (auto &particle : particle_handler)
    {
      auto &particle_properties = particle.get_properties();
#if DEAL_II_VERSION_GTE(10, 0, 0)
      MOI[particle.get_local_index()] =
#else
      MOI[particle.get_id()] =
#endif
        0.1 * particle_properties[DEM::PropertiesIndex::mass] *
        particle_properties[DEM::PropertiesIndex::dp] *
        particle_properties[DEM::PropertiesIndex::dp];
    }
}

template <int dim>
std::shared_ptr<Integrator<dim>>
CFDDEMSolver<dim>::set_integrator_type()
{
  if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
        .integration_method == Parameters::Lagrangian::ModelParameters::
                                 IntegrationMethod::velocity_verlet)
    {
      integrator_object = std::make_shared<VelocityVerletIntegrator<dim>>();
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .integration_method == Parameters::Lagrangian::ModelParameters::
                                      IntegrationMethod::explicit_euler)
    {
      integrator_object = std::make_shared<ExplicitEulerIntegrator<dim>>();
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .integration_method ==
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
std::shared_ptr<PPContactForce<dim>>
CFDDEMSolver<dim>::set_pp_contact_force()
{
  if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
        .pp_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PPContactForceModel::pp_linear)
    {
      pp_contact_force_object =
        std::make_shared<PPLinearForce<dim>>(dem_parameters);
    }
  else if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
             .pp_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::PPContactForceModel::
             pp_nonlinear)
    {
      pp_contact_force_object =
        std::make_shared<PPNonLinearForce<dim>>(dem_parameters);
    }
  else
    {
      throw "The chosen particle-particle contact force model is invalid";
    }
  return pp_contact_force_object;
}

template <int dim>
std::shared_ptr<PWContactForce<dim>>
CFDDEMSolver<dim>::set_pw_contact_force()
{
  std::vector<types::boundary_id> boundary_index =
    this->triangulation->get_boundary_ids();

  if (this->cfd_dem_simulation_parameters.dem_parameters.model_parameters
        .pw_contact_force_method ==
      Parameters::Lagrangian::ModelParameters::PWContactForceModel::pw_linear)
    {
      pw_contact_force_object = std::make_shared<PWLinearForce<dim>>(
        dem_parameters.boundary_conditions.boundary_translational_velocity,
        dem_parameters.boundary_conditions.boundary_rotational_speed,
        dem_parameters.boundary_conditions.boundary_rotational_vector,
        triangulation_cell_diameter,
        dem_parameters,
        boundary_index);
    }
  else if (dem_parameters.model_parameters.pw_contact_force_method ==
           Parameters::Lagrangian::ModelParameters::PWContactForceModel::
             pw_nonlinear)
    {
      pw_contact_force_object = std::make_shared<PWNonLinearForce<dim>>(
        dem_parameters.boundary_conditions.boundary_translational_velocity,
        dem_parameters.boundary_conditions.boundary_rotational_speed,
        dem_parameters.boundary_conditions.boundary_rotational_vector,
        triangulation_cell_diameter,
        dem_parameters,
        boundary_index);
    }
  else
    {
      throw "The chosen particle-wall contact force model is invalid";
    }
  return pw_contact_force_object;
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

#if DEAL_II_VERSION_GTE(10, 0, 0)
      types::particle_index particle_id = particle->get_local_index();
#else
      types::particle_index particle_id = particle->get_id();
#endif

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
CFDDEMSolver<dim>::dem_iterator(unsigned int counter)
{
  // dem_contact_build carries out the particle-particle and particle-wall
  // broad and fine searches, sort_particles_into_subdomains_and_cells, and
  // exchange_ghost
  dem_contact_build(counter);

  // Particle-particle contact force
  pp_contact_force_object->calculate_pp_contact_force(local_adjacent_particles,
                                                      ghost_adjacent_particles,
                                                      dem_time_step,
                                                      momentum,
                                                      force);

  // Particles-walls contact force:
  particle_wall_contact_force();

  // Add fluid-particle interaction force to the force container
  add_fluid_particle_interaction_force();

  // Integration correction step (after force calculation)
  // In the first step, we have to obtain location of particles at half-step
  // time
  if (this->simulation_control->get_step_number() == 0)
    {
      integrator_object->integrate_half_step_location(
        this->particle_handler,
        dem_parameters.lagrangian_physical_properties.g,
        force,
        dem_time_step);
    }
  else
    {
      integrator_object->integrate(
        this->particle_handler,
        dem_parameters.lagrangian_physical_properties.g,
        dem_time_step,
        momentum,
        force,
        MOI);
    }
}

template <int dim>
void
CFDDEMSolver<dim>::dem_contact_build(unsigned int counter)
{
  // Check to see if it is contact search step
  contact_detection_step = (this->*check_contact_search_step)(counter);

  // Sort particles in cells
  if (contact_detection_step || checkpoint_step || load_balance_step)
    {
      this->pcout << "DEM contact search at dem step " << counter << std::endl;

      this->particle_handler.sort_particles_into_subdomains_and_cells();

#if DEAL_II_VERSION_GTE(10, 0, 0)
      displacement.resize(
        this->particle_handler.get_max_local_particle_index());
#else
      {
        unsigned int max_particle_id = 0;
        for (const auto &particle : this->particle_handler)
          max_particle_id = std::max(max_particle_id, particle.get_id());
        displacement.resize(max_particle_id + 1);
      }
#endif
      force.resize(displacement.size());
      momentum.resize(displacement.size());

      this->particle_handler.exchange_ghost_particles(true);

      // Updating moment of inertia container
      update_moment_of_inertia(this->particle_handler, MOI);
    }
  else
    {
      this->particle_handler.update_ghost_particles();
    }

  // Broad particle-particle contact search
  // TODO add checkpoint step
  if (load_balance_step || checkpoint_step || contact_detection_step)
    {
      pp_broad_search_object.find_particle_particle_contact_pairs(
        this->particle_handler,
        &cells_local_neighbor_list,
        &cells_ghost_neighbor_list,
        local_contact_pair_candidates,
        ghost_contact_pair_candidates);


      // Particle-wall broad contact search
      particle_wall_broad_search();

      localize_contacts<dim>(&local_adjacent_particles,
                             &ghost_adjacent_particles,
                             &pw_pairs_in_contact,
                             &pfw_pairs_in_contact,
                             local_contact_pair_candidates,
                             ghost_contact_pair_candidates,
                             pw_contact_candidates,
                             pfw_contact_candidates);


      locate_local_particles_in_cells<dim>(this->particle_handler,
                                           particle_container,
                                           ghost_adjacent_particles,
                                           local_adjacent_particles,
                                           pw_pairs_in_contact,
                                           pfw_pairs_in_contact,
                                           particle_points_in_contact,
                                           particle_lines_in_contact);

      // Particle-particle fine search
      pp_fine_search_object.particle_particle_fine_search(
        local_contact_pair_candidates,
        ghost_contact_pair_candidates,
        local_adjacent_particles,
        ghost_adjacent_particles,
        particle_container,
        neighborhood_threshold_squared);

      // Particles-wall fine search
      particle_wall_fine_search();

      // Reset different steps
      checkpoint_step   = false;
      load_balance_step = false;
    }

  // Visualization
  if (this->simulation_control->is_output_iteration())
    {
      write_DEM_output_results();
    }

  // TODO add DEM post-processing

  // TODO checkpointing should be defined for CFD-DEM. We have to write
  // checkpointing for the GLS VANS solver and use the DEM checkpoint that
  // already exists.
}

template <int dim>
void
CFDDEMSolver<dim>::write_DEM_output_results()
{
  const std::string folder = this->cfd_dem_simulation_parameters.dem_parameters
                               .simulation_control.output_folder;
  const std::string particles_solution_name =
    this->cfd_dem_simulation_parameters.dem_parameters.simulation_control
      .output_name +
    "particles";
  const unsigned int iter = this->simulation_control->get_step_number();
  const double       time = this->simulation_control->get_current_time();
  const unsigned int group_files =
    this->cfd_dem_simulation_parameters.dem_parameters.simulation_control
      .group_files;

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
CFDDEMSolver<dim>::particle_wall_broad_search()
{
  // Particle - wall contact candidates
  pw_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cell_object.get_boundary_cells_information(),
    this->particle_handler,
    pw_contact_candidates);

  // Particle - floating wall contact pairs
  if (this->cfd_dem_simulation_parameters.dem_parameters.floating_walls
        .floating_walls_number > 0)
    {
      pw_broad_search_object.find_particle_floating_wall_contact_pairs(
        boundary_cell_object.get_boundary_cells_with_floating_walls(),
        this->particle_handler,
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        pfw_contact_candidates);
    }

  particle_point_contact_candidates =
    particle_point_line_broad_search_object.find_particle_point_contact_pairs(
      this->particle_handler,
      boundary_cell_object.get_boundary_cells_with_points());

  if (dim == 3)
    {
      particle_line_contact_candidates =
        particle_point_line_broad_search_object
          .find_particle_line_contact_pairs(
            this->particle_handler,
            boundary_cell_object.get_boundary_cells_with_lines());
    }
}

template <int dim>
void
CFDDEMSolver<dim>::particle_wall_fine_search()
{
  // Particle - wall fine search
  pw_fine_search_object.particle_wall_fine_search(pw_contact_candidates,
                                                  pw_pairs_in_contact);

  // Particle - floating wall fine search
  if (this->cfd_dem_simulation_parameters.dem_parameters.floating_walls
        .floating_walls_number > 0)
    {
      pw_fine_search_object.particle_floating_wall_fine_search(
        pfw_contact_candidates,
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        pfw_pairs_in_contact);
    }

  particle_points_in_contact =
    particle_point_line_fine_search_object.particle_point_fine_search(
      particle_point_contact_candidates, neighborhood_threshold_squared);

  if (dim == 3)
    {
      particle_lines_in_contact =
        particle_point_line_fine_search_object.particle_line_fine_search(
          particle_line_contact_candidates, neighborhood_threshold_squared);
    }
}

template <int dim>
void
CFDDEMSolver<dim>::particle_wall_contact_force()
{
  // Particle-wall contact force
  pw_contact_force_object->calculate_pw_contact_force(pw_pairs_in_contact,
                                                      dem_time_step,
                                                      momentum,
                                                      force);

  if (this->cfd_dem_simulation_parameters.dem_parameters.forces_torques
        .calculate_force_torque)
    {
      forces_boundary_information[this->simulation_control->get_step_number()] =
        pw_contact_force_object->get_force();
      torques_boundary_information[this->simulation_control
                                     ->get_step_number()] =
        pw_contact_force_object->get_torque();
    }

  // Particle-floating wall contact force
  if (this->cfd_dem_simulation_parameters.dem_parameters.floating_walls
        .floating_walls_number > 0)
    {
      pw_contact_force_object->calculate_pw_contact_force(pfw_pairs_in_contact,
                                                          dem_time_step,
                                                          momentum,
                                                          force);
    }

  particle_point_line_contact_force_object
    .calculate_particle_point_contact_force(
      &particle_points_in_contact,
      dem_parameters.lagrangian_physical_properties,
      force);

  if (dim == 3)
    {
      particle_point_line_contact_force_object
        .calculate_particle_line_contact_force(
          &particle_lines_in_contact,
          dem_parameters.lagrangian_physical_properties,
          force);
    }
}


template <int dim>
void
CFDDEMSolver<dim>::solve()
{
  // This is enforced to 1 right now because it does not provide
  // better speed-up than using MPI. This could be eventually changed...
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  // Reading DEM start file information
  if (this->cfd_dem_simulation_parameters.void_fraction->read_dem == true &&
      this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
          .restart == false)
    read_dem();

  this->setup_dofs();
  this->calculate_void_fraction(this->simulation_control->get_current_time());
  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  // Initilize DEM parameters
  initialize_dem_parameters();

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_control->is_at_start())
        {
          this->first_iteration();

          if (this->cfd_dem_simulation_parameters.cfd_parameters.test.enabled)
            { // Write particle Velocity
              for (auto &particle : this->particle_handler)
                {
                  auto particle_properties = particle.get_properties();
                  this->pcout
                    << "Particle Summary"
                    << "\n"
                    << "--------------------------------------------------------------------------"
                    << "--------------------------------------------------------------------------"
                    << "\n"
                    << "id: " << particle.get_id() << ",  "
                    << "x: " << particle.get_location()[0] << ",  "
                    << "y: " << particle.get_location()[1] << ",  "
                    << "z: " << particle.get_location()[2] << ",  "
                    << "v_x: " << particle_properties[DEM::PropertiesIndex::v_x]
                    << ",  "
                    << "v_y: " << particle_properties[DEM::PropertiesIndex::v_y]
                    << ",  "
                    << "vz: " << particle_properties[DEM::PropertiesIndex::v_z]
                    << "\n"
                    << "--------------------------------------------------------------------------"
                    << "--------------------------------------------------------------------------"
                    << std::endl;
                }
            }

          this->pcout << "Starting DEM iterations at step "
                      << this->simulation_control->get_step_number()
                      << std::endl;
          for (unsigned int dem_counter = 0; dem_counter < coupling_frequency;
               ++dem_counter)
            {
              TimerOutput::Scope t(this->computing_timer, "DEM_Iterator");
              // dem_iterator carries out the particle-particle and
              // particle_wall force calculations, integration and
              // update_ghost
              dem_iterator(dem_counter);
            }

          this->pcout << "Finished " << coupling_frequency << " DEM iterations "
                      << std::endl;
        }
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->iterate();

          if (this->cfd_dem_simulation_parameters.cfd_parameters.test.enabled)
            {
              // Write particle Velocity
              for (auto &particle : this->particle_handler)
                {
                  auto particle_properties = particle.get_properties();
                  this->pcout
                    << "Particle Summary"
                    << "\n"
                    << "--------------------------------------------------------------------------"
                    << "--------------------------------------------------------------------------"
                    << "\n"
                    << "id: " << particle.get_id() << ",  "
                    << "x: " << particle.get_location()[0] << ",  "
                    << "y: " << particle.get_location()[1] << ",  "
                    << "z: " << particle.get_location()[2] << ",  "
                    << "v_x: " << particle_properties[DEM::PropertiesIndex::v_x]
                    << ",  "
                    << "v_y: " << particle_properties[DEM::PropertiesIndex::v_y]
                    << ",  "
                    << "vz: " << particle_properties[DEM::PropertiesIndex::v_z]
                    << "\n"
                    << "--------------------------------------------------------------------------"
                    << "--------------------------------------------------------------------------"
                    << std::endl;
                }
            }

          this->pcout << "Starting DEM iterations at step "
                      << this->simulation_control->get_step_number()
                      << std::endl;

          for (unsigned int dem_counter = 0; dem_counter < coupling_frequency;
               ++dem_counter)
            {
              TimerOutput::Scope t(this->computing_timer, "DEM_Iterator");
              // dem_iterator carries out the particle-particle and
              // particle_wall force calculations, integration and
              // update_ghost
              dem_iterator(dem_counter);
            }
          this->pcout << "Finished " << coupling_frequency << " DEM iterations "
                      << std::endl;
        }

      this->postprocess(false);
      this->finish_time_step();

      if (this->cfd_dem_simulation_parameters.cfd_dem.post_processing)
        this->post_processing();

      // Load balancing
      load_balance_step = (this->*check_load_balance_step)(0);
    }
  this->finish_simulation();
}

// Pre-compile the 2D and 3D CFD-DEM solver to ensure that the
// library is valid before we actually compile the solver This greatly
// helps with debugging
template class CFDDEMSolver<2>;
template class CFDDEMSolver<3>;
