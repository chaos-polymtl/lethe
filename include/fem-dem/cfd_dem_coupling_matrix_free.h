// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cfd_dem_coupling_h
#define lethe_cfd_dem_coupling_h

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/dem_action_manager.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/periodic_boundaries_manipulator.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/fluid_dynamics_vans_matrix_free.h>


using namespace dealii;

/**
 * @brief Solver class for coupled Computational Fluid Dynamics and Discrete Element Method (CFD-DEM) simulations using a matrix-free formulation.
 *
 * This class implements a two-way coupled CFD-DEM solver that combines the
 * Volume-averaged Navier-Stokes (VANS) equations for fluid flow with the
 * soft-sphere Discrete Element Method for particle dynamics. The solver
 * enables simulation of fluid-particle interactions in systems containing
 * dense particle suspensions.
 *
 * The CFD component uses Galerkin Least Squares (GLS) stabilization for
 * the VANS equations, while the DEM component handles particle-particle
 * and particle-wall collisions using the soft-sphere contact model. A
 * matrix-free formulation is used to solve the VANS equation.
 *
 * @tparam dim Spatial dimension of the simulation (2 or 3).
 * While 2D simulation are theoretically supported, they have not been tested.
 *
 * @ingroup solvers
 */
template <int dim>
class CFDDEMMatrixFree : public FluidDynamicsVANSMatrixFree<dim>
{
public:
  /**
   * @brief Constructor for the CFD-DEM matrix-free solver.
   *
   * @param[in] param Simulation parameters containing all CFD-DEM
   * configuration options
   */
  CFDDEMMatrixFree(CFDDEMSimulationParameters<dim> &param);



  /**
   * @brief Engine of the CFD-DEM solver. Calls all the necessary functions to
   * set parameters, solve the simulation, and finish the simulation.
   */
  void
  solve() override;

private:
  /**
   * @brief Set up the various parameters related to the DEM simulation.
   *
   * Initializes DEM-specific parameters including contact detection methods,
   * integration schemes, and material properties based on the simulation
   * parameters provided during construction.
   */
  void
  dem_setup_parameters();

  /**
   * @brief Initialize the particle size distribution and related parameters.
   *
   * Sets up the particle size distribution type for the simulation and
   * calculates the maximum particle diameter and neighborhood threshold
   * squared values required for contact detection algorithms.
   */
  void
  setup_distribution_type();

  /**
   * @brief Set up the time integration method for particle motion.
   *
   * Creates and configures the appropriate integration scheme for solving
   * particle equations of motion (e.g., velocity Verlet, explicit Euler).
   *
   * @return Shared pointer to the configured integration object
   */
  std::shared_ptr<Integrator<dim, DEM::CFDDEMProperties::PropertiesIndex>>
  set_integrator_type();

  /**
   * @brief Initialize fundamental DEM simulation parameters.
   *
   * Sets up basic DEM parameters including time stepping, contact detection
   * frequency, and other core simulation settings.
   */
  void
  initialize_dem_parameters();

  /**
   * @brief Read DEM restart files from a previous simulation.
   *
   * Loads particle data, contact information, and other DEM state variables
   * from checkpoint files to continue a previously interrupted simulation.
   */
  void
  read_dem() override;

  /**
   * @brief Set the initial condition. The version in this class does not output
   * the solution, since the void fraction is not initialized yet.
   *
   * @param[in] initial_condition_type Type of method  use to impose initial
   *condition.
   *
   * @param[in] restart Indicator if the simulation is being restarted or not.
   *
   **/

  void
  set_initial_condition(
    const Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    const bool restart = false) override
  {
    unsigned int ref_iter = 0;
    do
      {
        if (ref_iter > 0)
          this->refine_mesh();

        this->set_initial_condition_fd(initial_condition_type, restart);
        if (!restart)
          {
            this->multiphysics->set_initial_conditions();
          }
        ref_iter++;
      }
    while (
      ref_iter <
        (this->simulation_parameters.mesh_adaptation.initial_refinement + 1) &&
      restart == false);

    this->pcout
      << "---------------------------------------------------------------"
      << std::endl;
  }

  /**
   * @brief Gather TableHandler objects for serialization during checkpointing.
   *
   * Collects all TableHandler objects that need to be serialized and
   * deserialized during checkpoint operations, including their associated
   * file names for proper storage and retrieval.
   *
   * @return Vector of OutputStructTableHandler containing references to
   * TableHandler objects and their corresponding file names
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Write CFD-DEM simulation checkpoint files.
   *
   * Saves the complete state of the CFD-DEM simulation including fluid
   * variables, particle data, and auxiliary information to enable
   * simulation restart from the current time step.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Read CFD-DEM simulation checkpoint files.
   *
   * Loads the complete simulation state from previously written checkpoint
   * files to restart a CFD-DEM simulation from a specific time step.
   */
  void
  read_checkpoint() override;

  /**
   * @brief Check if contact detection is required and execute it if necessary.
   *
   * Determines whether a contact detection step should be performed based on
   * the contact detection frequency and executes the appropriate contact
   * detection algorithm.
   *
   * @param[in] counter Current DEM iteration number
   */
  void
  check_contact_detection_method(unsigned int counter);

  /**
   * @brief Check if load balancing is required and perform it if necessary.
   *
   * Evaluates the current load distribution across MPI processes and
   * performs load balancing to maintain computational efficiency when
   * particle distribution becomes uneven.
   */
  void
  load_balance();

  /**
   * @brief Check if particles need to be inserted and perform insertion.
   *
   * Determines if new particles should be inserted into the simulation
   * domain based on the insertion parameters and executes particle
   * insertion using the configured insertion method.
   */
  void
  insert_particles();


  /**
   * @brief Calculate and add fluid-particle interaction forces and torques.
   *
   * Computes the hydrodynamic forces acting on particles due to fluid motion
   * and adds them to the particle force container for use in DEM integration.
   */
  void
  add_fluid_particle_interaction();

  /**
   * @brief Calculate contact forces between particles and walls.
   *
   * Computes collision forces and torques resulting from particle-wall
   * interactions using the soft-sphere contact model.
   */
  void
  particle_wall_contact_force();

  /**
   * @brief Write DEM simulation results to output files.
   *
   * Generates visualization output files containing particle positions,
   * velocities, forces, and other relevant DEM variables in parallel
   * VTU format for post-processing.
   */
  void
  write_dem_output_results();

  /**
   * @brief Calculate and report particle statistics to the terminal.
   *
   * Computes and displays comprehensive statistics about particle behavior
   * including minimum, maximum, and average values for linear velocity,
   * angular velocity, kinetic energies, and contact search performance
   * metrics.
   */
  void
  report_particle_statistics();

  /**
   * @brief Perform fluid dynamics post-processing after each iteration.
   *
   * Executes post-processing operations specific to the fluid dynamics
   * component including output writing and monitoring.
   *
   * @param[in] first_iteration Flag indicating if this is the first iteration
   */
  void
  postprocess_fd(bool first_iteration) override;

  /**
   * @brief Perform CFD-DEM coupled post-processing after each iteration.
   *
   * Executes post-processing operations specific to the coupled CFD-DEM
   * simulation including coupled field analysis and specialized output
   * generation.
   */
  void
  postprocess_cfd_dem();

  /**
   * @brief Calculate dynamic flow control accounting for void fraction.
   *
   * Performs dynamic flow control calculations that incorporate the local
   * void fraction (porosity) field in the average velocity computation,
   * enabling proper flow control in particle-laden flows.
   */
  void
  dynamic_flow_control() override;

  /**
   * @brief Sort particles into subdomains and cells for efficient processing.
   *
   * Redistributes particles among computational subdomains and associates
   * them with their containing cells, then reinitialize containers that
   * depend on local particle identifiers.
   */
  void
  sort_particles_into_subdomains_and_cells();

  /**
   * @brief Execute a DEM time step.
   *
   * Performs a full DEM time step including particle-particle and
   * particle-wall contact force calculations, time integration of
   * particle motion, and updates to particle positions and velocities.
   *
   * @param[in] counter Current DEM iteration number
   */
  void
  dem_iterator(unsigned int counter);

  /**
   * @brief Calculate the new DEM time-step.
   *
   * Calculate the new DEM time-step in a CFD-DEM simulation with adaptive
   * time-stepping
   *
   */
  void
  update_dem_time_step()
  {
    dem_time_step =
      this->simulation_control->get_time_step() / coupling_frequency;
    const double time_step_rayleigh_ratio = dem_time_step / rayleigh_time_step;
    this->pcout << "DEM time-step is " << time_step_rayleigh_ratio * 100
                << "% of Rayleigh time step" << std::endl;
  }

  /**
   * @brief Perform contact detection and build contact lists.
   *
   * Checks if a contact search should be performed based on the contact
   * detection frequency and executes the appropriate contact detection
   * algorithm to identify particle-particle and particle-wall contacts.
   *
   * @param[in] counter Current DEM iteration number
   */
  void
  dem_contact_build(unsigned int counter);

  /// Frequency of coupling between CFD and DEM solvers (in DEM time steps)
  unsigned int coupling_frequency;

  /// Gravitational acceleration vector
  // TODO - Refactor the gravity to not be a dangling Tensor in the class (BB)
  Tensor<1, 3> g;

  /// Container for particle interaction outcomes (forces and torques)
  ParticleInteractionOutcomes<DEM::CFDDEMProperties::PropertiesIndex>
    contact_outcome;

  /// Displacement vector for particle motion calculations
  std::vector<double> displacement;

  /// Moment of inertia values for all particles
  std::vector<double> MOI;

  /// Squared threshold distance for neighborhood particle detection
  double neighborhood_threshold_squared;

  /// Container for particle size distribution objects
  std::vector<std::shared_ptr<Distribution>> size_distribution_object_container;

  /// Maximum diameter among all particles in the simulation
  double maximum_particle_diameter;

  /// Smallest criterion used for contact search algorithms
  double smallest_contact_search_criterion;

  /**
   * @brief Object to store ongoing collision information.
   *
   * Stores particle ID, boundary ID, velocity tensor, angular velocity
   * tensor, and collision time for collisions currently in progress.
   */
  OngoingCollisionLog<dim> ongoing_collision_log;

  /**
   * @brief Object to store completed collision events.
   *
   * Stores complete collision event data including start and end information
   * for post-processing and analysis purposes.
   */
  CompletedCollisionLog<dim> collision_event_log;

  /// Manager for particle-particle and particle-wall contact detection and
  /// handling
  DEMContactManager<dim, DEM::CFDDEMProperties::PropertiesIndex>
    contact_manager;

  /// Load balancing handler for distributing particles across MPI processes
  LagrangianLoadBalancing<dim, DEM::CFDDEMProperties::PropertiesIndex>
    load_balancing;

  /// Handler for point and line contact forces between particles
  ParticlePointLineForce<dim, DEM::CFDDEMProperties::PropertiesIndex>
    particle_point_line_contact_force_object;

  /// Time integration scheme for particle motion equations
  std::shared_ptr<Integrator<dim, DEM::CFDDEMProperties::PropertiesIndex>>
    integrator_object;

  /// Particle insertion method handler
  std::shared_ptr<Insertion<dim, DEM::CFDDEMProperties::PropertiesIndex>>
    insertion_object;

  /// Contact force model for particle-particle interactions
  std::shared_ptr<
    ParticleParticleContactForceBase<dim,
                                     DEM::CFDDEMProperties::PropertiesIndex>>
    particle_particle_contact_force_object;

  /// Contact force model for particle-wall interactions
  std::shared_ptr<
    ParticleWallContactForceBase<dim, DEM::CFDDEMProperties::PropertiesIndex>>
    particle_wall_contact_force_object;

  /// Visualization handler for DEM output generation
  Visualization<dim, DEM::CFDDEMProperties::PropertiesIndex>
    visualization_object;

  /// Information about boundary cells for contact detection
  BoundaryCellsInformation<dim> boundary_cell_object;

  /// Mesh information for solid surfaces in contact with particles
  typename dem_data_structures<dim>::solid_surfaces_mesh_information
    solid_surfaces_mesh_info;

  /// Updated boundary points and their corresponding normal vectors
  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;

  /// Force information stored on boundary elements
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;

  /// Torque information stored on boundary elements
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;

  /// Information about cells involved in periodic boundary conditions
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information;

  /// Handler for periodic boundary condition manipulations
  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;

  /// Object for handling adaptive sparse contact detection algorithms
  AdaptiveSparseContacts<dim, DEM::CFDDEMProperties::PropertiesIndex>
    sparse_contacts_object;

  /// Counter for contact searches performed in the current CFD iteration
  unsigned int contact_search_counter;

  /// Total number of contact searches since simulation start
  unsigned int contact_search_total_number;

  /// Statistics container for timing and contact list information
  statistics contact_list;

  /// Particle properties handler for DEM simulation
  DEM::ParticleProperties<dim, DEM::CFDDEMProperties::PropertiesIndex>
    properties_class;

  /// PVD handler for grid visualization output
  PVDHandler grid_pvdhandler;

  /// PVD handler for particle visualization output
  PVDHandler particles_pvdhandler;

  /// DEM simulation parameters container
  DEMSolverParameters<dim> dem_parameters;

  /// Time step size for DEM integration
  double dem_time_step;

  /// Rayleigh time step
  double rayleigh_time_step;

  /// TableHandler for post-processing phase volume calculations
  /// Stores total fluid volume and total particle volume data
  TableHandler table_phase_volumes;

  /// Pointer to DEM action manager for handling DEM-specific operations
  DEMActionManager *dem_action_manager;
};
#endif
