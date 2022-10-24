#include <core/solutions_output.h>

#include <dem/dem_solver_parameters.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/find_maximum_particle_size.h>
#include <dem/gear3_integrator.h>
#include <dem/particle_particle_linear_force.h>
#include <dem/particle_particle_nonlinear_force.h>
#include <dem/particle_wall_linear_force.h>
#include <dem/particle_wall_nonlinear_force.h>
#include <dem/post_processing.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <fem-dem/cfd_dem_coupling.h>

#include <fstream>

template <int dim>
bool
check_contact_detection_method(
  unsigned int                          counter,
  CFDDEMSimulationParameters<dim> &     param,
  std::vector<double> &                 displacement,
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
    { // The sorting into subdomain step checks whether or not the current time
      // step
      // is a step that requires sorting particles into subdomains and cells.
      // This is applicable if any of the following three conditions apply:if
      // its a load balancing step, a restart simulation step, or a contact
      // detection tsep.
      bool sorting_in_subdomains_step =
        (checkpoint_step || load_balance_step || contact_detection_step);

      if (sorting_in_subdomains_step)
        displacement.resize(particle_handler.get_max_local_particle_index());

      contact_detection_step =
        find_contact_detection_step<dim>(particle_handler,
                                         simulation_control->get_time_step() /
                                           param.cfd_dem.coupling_frequency,
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
  CFDDEMSimulationParameters<dim> &     param,
  Particles::ParticleHandler<dim, dim> &particle_handler,
  const MPI_Comm &                      mpi_communicator,
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
  , this_mpi_process(Utilities::MPI::this_mpi_process(this->mpi_communicator))
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(this->mpi_communicator))
{
  coupling_frequency =
    this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency;

  standard_deviation_multiplier = 2.5;

  // Initialize DEM Parameters
  dem_parameters.lagrangian_physical_properties =
    this->cfd_dem_simulation_parameters.dem_parameters
      .lagrangian_physical_properties;
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


  maximum_particle_diameter =
    find_maximum_particle_size(dem_parameters.lagrangian_physical_properties,
                               standard_deviation_multiplier);
  neighborhood_threshold_squared =
    std::pow(dem_parameters.model_parameters.neighborhood_threshold *
               maximum_particle_diameter,
             2);

  if (dem_parameters.mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        tria,
        this->cfd_dem_simulation_parameters.dem_parameters.mesh.grid_type,
        this->cfd_dem_simulation_parameters.dem_parameters.mesh.grid_arguments);
    }
  else if (dem_parameters.mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(tria);
      std::ifstream input_file(dem_parameters.mesh.file_name);
      grid_in.read_msh(input_file);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  triangulation_cell_diameter = 0.5 * GridTools::diameter(tria);

  //   Finding the smallest contact search frequency criterion between (smallest
  //   cell size - largest particle radius) and (security factor * (blab
  //   diamater - 1) *  largest particle radius). This value is used in
  //   find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(tria) -
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
    }

  // Initilize contact detection step
  contact_detection_step = false;
  load_balance_step      = false;

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
  TimerOutput::Scope timer(this->computing_timer, "Write_Checkpoint");
  this->pcout << "Writing restart file" << std::endl;

  std::string prefix = this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);
      this->pvdhandler.save(prefix);
      particles_pvdhandler.save(prefix_particles);

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

  // Prepare particle handler for serialization
  this->particle_handler.prepare_for_serialization();

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

  if (this->simulation_parameters.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);

      this->flow_rate = calculate_flow_rate(
        this->dof_handler,
        this->present_solution,
        this->simulation_parameters.flow_control.boundary_flow_id,
        *this->face_quadrature,
        *this->mapping);
    }

  this->multiphysics->read_checkpoint();

  // Deserialize particles have the triangulation has been read
  this->particle_handler.deserialize();
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

  // Prepare particle handle for serialization
  this->particle_handler.prepare_for_coarsening_and_refinement();

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

  if (dem_parameters.boundary_conditions.BC_type ==
      Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
    {
      periodic_boundaries_object.map_periodic_cells(*parallel_triangulation);
    }

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

  // Unpack particle handler after load balancing step
  this->particle_handler.unpack_after_coarsening_and_refinement();
}

template <int dim>
void
CFDDEMSolver<dim>::initialize_dem_parameters()
{
  this->pcout << "Initializing DEM parameters " << std::endl;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  // Finding cell neighbors
  cell_neighbors_object.find_cell_neighbors(*parallel_triangulation,
                                            cells_local_neighbor_list,
                                            cells_ghost_neighbor_list);

  if (dem_parameters.boundary_conditions.BC_type ==
      Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
    {
      periodic_boundaries_object.set_periodic_boundaries_information(
        dem_parameters.boundary_conditions.outlet_boundaries,
        dem_parameters.boundary_conditions.periodic_boundaries,
        dem_parameters.boundary_conditions.periodic_direction);

      periodic_boundaries_object.map_periodic_cells(*parallel_triangulation);
    }

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
    *parallel_triangulation,
    triangulation_cell_diameter);

  this->particle_handler.sort_particles_into_subdomains_and_cells();

  displacement.resize(this->particle_handler.get_max_local_particle_index());

  force.resize(displacement.size());
  torque.resize(displacement.size());


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
  MOI.resize(torque.size());

  for (auto &particle : particle_handler)
    {
      auto &particle_properties = particle.get_properties();
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
CFDDEMSolver<dim>::dem_iterator(unsigned int counter)
{
  // dem_contact_build carries out the particle-particle and particle-wall
  // broad and fine searches, sort_particles_into_subdomains_and_cells, and
  // exchange_ghost
  dem_contact_build(counter);

  // Particle-particle contact force
  particle_particle_contact_force_object
    ->calculate_particle_particle_contact_force(local_adjacent_particles,
                                                ghost_adjacent_particles,
                                                dem_time_step,
                                                torque,
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
        dem_time_step,
        torque,
        force,
        MOI);
    }
  else
    {
      integrator_object->integrate(
        this->particle_handler,
        dem_parameters.lagrangian_physical_properties.g,
        dem_time_step,
        torque,
        force,
        MOI);
    }

  // Particles displacement if passing through a periodic boundary
  periodic_boundaries_object.execute_particles_displacement(
    this->particle_handler);
}

template <int dim>
void
CFDDEMSolver<dim>::dem_contact_build(unsigned int counter)
{
  // Check to see if it is contact search step
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
  // The last statement in the if condition allows the calling of
  // dem_contact_build for the first DEM iteration in the first CFD time step
  // directly after reading the dem initial checkpoint files

  if (contact_detection_step || checkpoint_step || load_balance_step ||
      (this->simulation_control->is_at_start() && (counter == 0)))
    {
      this->pcout << "DEM contact search at dem step " << counter << std::endl;

      this->particle_handler.sort_particles_into_subdomains_and_cells();

      displacement.resize(
        this->particle_handler.get_max_local_particle_index());
      force.resize(displacement.size());
      torque.resize(displacement.size());

      this->particle_handler.exchange_ghost_particles(true);

      // Updating moment of inertia container
      update_moment_of_inertia(this->particle_handler, MOI);
    }
  else
    {
      this->particle_handler.update_ghost_particles();
    }

  // Broad particle-particle contact search
  if (load_balance_step || checkpoint_step || contact_detection_step ||
      (this->simulation_control->is_at_start() && (counter == 0)))
    {
      particle_particle_broad_search_object
        .find_particle_particle_contact_pairs(this->particle_handler,
                                              cells_local_neighbor_list,
                                              cells_ghost_neighbor_list,
                                              local_contact_pair_candidates,
                                              ghost_contact_pair_candidates);


      // Particle-wall broad contact search
      particle_wall_broad_search();

      localize_contacts<dim>(local_adjacent_particles,
                             ghost_adjacent_particles,
                             particle_wall_in_contact,
                             particle_floating_wall_in_contact,
                             particle_floating_mesh_in_contact,
                             local_contact_pair_candidates,
                             ghost_contact_pair_candidates,
                             particle_wall_candidates,
                             particle_floating_wall_candidates,
                             particle_floating_mesh_candidates);


      locate_local_particles_in_cells<dim>(this->particle_handler,
                                           particle_container,
                                           ghost_adjacent_particles,
                                           local_adjacent_particles,
                                           particle_wall_in_contact,
                                           particle_floating_wall_in_contact,
                                           particle_floating_mesh_in_contact,
                                           particle_points_in_contact,
                                           particle_lines_in_contact);

      // Particle-particle fine search
      particle_particle_fine_search_object.particle_particle_fine_search(
        local_contact_pair_candidates,
        ghost_contact_pair_candidates,
        local_adjacent_particles,
        ghost_adjacent_particles,
        particle_container,
        neighborhood_threshold_squared);

      // Particles-wall fine search
      particle_wall_fine_search();

      // Reset different steps. The contact build should be performed everytime
      // we restart the simulation or everytime load balancing is performed. At
      // the end of the restart step or the load balance step, and after all
      // necessary contact build, vector resizing, and solution transfers have
      // been performed, this functions are set to false. The checkpoint_step
      // remains false for the duration of the simulation while the load
      // balancing is reset everytime load balancing is called.
      checkpoint_step   = false;
      load_balance_step = false;
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
CFDDEMSolver<dim>::particle_wall_broad_search()
{
  // Particle - wall contact candidates
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cell_object.get_boundary_cells_information(),
    this->particle_handler,
    particle_wall_candidates);

  // Particle - floating wall contact pairs
  if (dem_parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_broad_search_object
        .find_particle_floating_wall_contact_pairs(
          boundary_cell_object.get_boundary_cells_with_floating_walls(),
          this->particle_handler,
          dem_parameters.floating_walls,
          this->simulation_control->get_current_time(),
          particle_floating_wall_candidates);
    }

  particle_point_candidates =
    particle_point_line_broad_search_object.find_particle_point_contact_pairs(
      this->particle_handler,
      boundary_cell_object.get_boundary_cells_with_points());

  if (dim == 3)
    {
      particle_line_candidates =
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
  particle_wall_fine_search_object.particle_wall_fine_search(
    particle_wall_candidates, particle_wall_in_contact);

  // Particle - floating wall fine search
  if (dem_parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_fine_search_object.particle_floating_wall_fine_search(
        particle_floating_wall_candidates,
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        particle_floating_wall_in_contact);
    }

  particle_points_in_contact =
    particle_point_line_fine_search_object.particle_point_fine_search(
      particle_point_candidates, neighborhood_threshold_squared);

  if (dim == 3)
    {
      particle_lines_in_contact =
        particle_point_line_fine_search_object.particle_line_fine_search(
          particle_line_candidates, neighborhood_threshold_squared);
    }
}

template <int dim>
void
CFDDEMSolver<dim>::particle_wall_contact_force()
{
  // Particle-wall contact force
  particle_wall_contact_force_object->calculate_particle_wall_contact_force(
    particle_wall_in_contact, dem_time_step, torque, force);

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
        particle_floating_wall_in_contact, dem_time_step, torque, force);
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
CFDDEMSolver<dim>::dem_post_process_results()
{
  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  if (dem_parameters.post_processing.calculate_particles_average_velocity)
    {
      dem_post_processing_object.calculate_average_particles_velocity(
        *parallel_triangulation, this->particle_handler);

      dem_post_processing_object.write_average_particles_velocity(
        *parallel_triangulation,
        grid_pvdhandler,
        dem_parameters,
        this->simulation_control->get_current_time(),
        this->simulation_control->get_step_number(),
        this->mpi_communicator);
    }
  if (dem_parameters.post_processing.calculate_granular_temperature)
    {
      dem_post_processing_object.calculate_average_granular_temperature(
        *parallel_triangulation, this->particle_handler);

      dem_post_processing_object.write_granular_temperature(
        *parallel_triangulation,
        grid_pvdhandler,
        dem_parameters,
        this->simulation_control->get_current_time(),
        this->simulation_control->get_step_number(),
        this->mpi_communicator);
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
CFDDEMSolver<dim>::post_processing()
{
  if (this->cfd_dem_simulation_parameters.cfd_dem.post_processing)
    {
      this->GLSVANSSolver<dim>::post_processing();
      double particle_total_kinetic_energy =
        DEM::calculate_granular_kinetic_energy(this->particle_handler,
                                                     this->mpi_communicator).total;

      this->pcout << "Total particles kinetic energy: "
                  << particle_total_kinetic_energy << std::endl;
    }

  if (dem_parameters.post_processing.Lagrangian_post_processing)
    {
      if (this->simulation_control->get_step_number() >=
            dem_parameters.post_processing.initial_step &&
          this->simulation_control->get_step_number() <=
            dem_parameters.post_processing.end_step)
        {
          if (this->simulation_control->get_step_number() %
                dem_parameters.post_processing.output_frequency ==
              0)
            dem_post_process_results();
        }
    }
}

template <int dim>
void
CFDDEMSolver<dim>::print_particles_summary()
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
        << "v_x: " << particle_properties[DEM::PropertiesIndex::v_x] << ",  "
        << "v_y: " << particle_properties[DEM::PropertiesIndex::v_y] << ",  "
        << "vz: " << particle_properties[DEM::PropertiesIndex::v_z] << "\n"
        << "--------------------------------------------------------------------------"
        << "--------------------------------------------------------------------------"
        << std::endl;
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

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  // Initilize DEM parameters
  initialize_dem_parameters();

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          this->update_boundary_conditions();
        }

      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          this->vertices_cell_mapping();
          this->initialize_void_fraction();
          this->iterate();
        }
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
            refine_mesh();
          this->vertices_cell_mapping();
          this->calculate_void_fraction(
            this->simulation_control->get_current_time(), load_balance_step);
          this->iterate();
        }

      if (this->cfd_dem_simulation_parameters.cfd_parameters.test.enabled)
        {
          print_particles_summary();
        }

      this->pcout << "Starting DEM iterations at step "
                  << this->simulation_control->get_step_number() << std::endl;
      {
        TimerOutput::Scope t(this->computing_timer, "DEM_Iterator");

        for (unsigned int dem_counter = 0; dem_counter < coupling_frequency;
             ++dem_counter)
          {
            // dem_iterator carries out the particle-particle and
            // particle_wall force calculations, integration and
            // update_ghost
            dem_iterator(dem_counter);
          }
      }

      this->pcout << "Finished " << coupling_frequency << " DEM iterations "
                  << std::endl;

      this->postprocess(false);
      this->finish_time_step_fd();

      post_processing();

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
