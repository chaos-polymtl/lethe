// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_h
#define lethe_fluid_dynamics_vans_h

#include <core/parameters.h>

#include <solvers/fluid_dynamics_matrix_based.h>
#include <solvers/postprocessing_cfd.h>

#include <dem/dem.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/particle_projector.h>
#include <fem-dem/vans_assemblers.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

using namespace dealii;

/**
 * @brief Solver class for Volume-Averaged Navier-Stokes (VANS) equations with GLS stabilization.
 *
 * This class implements a computational fluid dynamics solver for the VANS
 * equations, which are used to model fluid flow in porous media or multiphase
 * systems where the presence of solid particles affects the fluid motion. The
 * solver uses Galerkin Least Squares (GLS) stabilization techniques including
 * SUPG (Streamline Upwind Petrov-Galerkin) and PSPG (Pressure Stabilizing
 * Petrov-Galerkin) methods.
 *
 * The VANS equations account for the volume fraction of the fluid phase (void
 * fraction) and incorporate particle-fluid interactions through appropriate
 * source terms. This solver forms the foundation for CFD-DEM simulations where
 * particle dynamics are coupled with fluid flow.
 *
 * @tparam dim Spatial dimension of the simulation (2 or 3)
 *
 * @ingroup solvers
 */

template <int dim>
class FluidDynamicsVANS : public FluidDynamicsMatrixBased<dim>
{
public:
  /**
   * @brief Constructor for the VANS solver.
   *
   * @param[in] nsparam CFD-DEM simulation parameters containing all
   * configuration options for the coupled fluid-particle simulation
   */
  FluidDynamicsVANS(CFDDEMSimulationParameters<dim> &nsparam);

  /**
   * @brief Destructor for the VANS solver.
   */
  ~FluidDynamicsVANS();

  /**
   * @brief Main solver loop for the VANS equations.
   *
   * Executes the complete solution procedure for the VANS equations including
   * time stepping, matrix assembly, linear system solution, and convergence
   * checking.
   */
  virtual void
  solve() override;

private:
  /**
   * @brief Assemble the diagonal mass matrix for the VANS system.
   *
   * Constructs the diagonal mass matrix used in time-dependent VANS
   * calculations, incorporating void fraction effects.
   *
   * @param[out] mass_matrix Sparse matrix to store the assembled mass matrix
   */
  void
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  /**
   * @brief Update solution vectors and constraint objects.
   *
   * Updates the solution vectors and applies constraints after each
   * iteration or time step in the VANS solution process.
   */
  void
  update_solution_and_constraints();

  /**
   * @brief Read DEM data for particle-fluid coupling.
   *
   * Loads particle information from DEM simulation for use in computing
   * void fraction and particle-fluid interaction terms.
   */
  virtual void
  read_dem();

protected:
  /**
   * @brief Set up degrees of freedom and initialize void fraction field.
   *
   * Associates degrees of freedom to each vertex of the finite elements
   * and initializes the void fraction field that represents the local
   * volume fraction available to the fluid phase.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Set the initial condition. The version in this class does not output
   * the solution, since the void fraction is not initialized yet.
   *
   * @param initial_condition_type Type of method  use to impose initial condition.
   *
   * @param restart Indicator if the simulation is being restarted or not.
   *
   **/

  void
  set_initial_condition(
    const Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    const bool                                          restart = false) override
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
  }

  /**
   * @brief Execute a single VANS iteration.
   *
   * Performs one iteration of the VANS solution procedure including
   * assembly, solution of the linear system, and solution update.
   */
  virtual void
  iterate() override;

  /**
   * @brief Calculate void fraction field based on particle positions.
   *
   * Computes the local void fraction (porosity) at each point in the
   * computational domain based on the current particle configuration.
   *
   * @param[in] time Current simulation time for time-dependent calculations
   */
  void
  calculate_void_fraction(const double time);

  /**
   * @brief Create mapping between vertices and cells.
   *
   * Establishes the relationship between mesh vertices and the cells
   * that contain them, which is needed for efficient void fraction
   * and particle-fluid coupling calculations.
   */
  void
  vertices_cell_mapping();

  /**
   * @brief Monitor mass conservation in the VANS system.
   *
   * Calculates and reports mass conservation metrics for the VANS
   * solution, taking into account the void fraction field.
   */
  virtual void
  monitor_mass_conservation();

  /**
   * @brief Complete operations at the end of a time step.
   *
   * Performs final operations required at the completion of each time
   * step including post-processing and time advancement for fluid dynamics.
   */
  virtual void
  finish_time_step_fd();

  /**
   * @brief Assemble the global system matrix for the VANS equations.
   *
   * Constructs the global system matrix by assembling contributions from
   * all cells, including fluid dynamics terms modified by void fraction
   * and particle-fluid interaction terms.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the global right-hand side vector for the VANS equations.
   *
   * Constructs the global right-hand side vector including all source terms,
   * boundary conditions, and particle-fluid coupling contributions.
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Assemble the local system matrix for a given cell.
   *
   * Computes the local contributions to the global system matrix for a
   * single cell, including VANS-specific terms with void fraction effects
   * and GLS stabilization terms. This function is used by WorkStream
   * for thread-safe parallel assembly.
   *
   * @param[in] cell Iterator pointing to the current active cell
   * @param[in,out] scratch_data Scratch data containing finite element
   * information at quadrature points for efficient assembly
   * @param[out] copy_data Copy data structure for storing local assembly
   * results before copying to global structures
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief Assemble the local right-hand side vector for a given cell.
   *
   * Computes the local contributions to the global right-hand side vector
   * for a single cell, including source terms, boundary conditions, and
   * particle-fluid interaction terms modified by void fraction.
   *
   * @param[in] cell Iterator pointing to the current active cell
   * @param[in,out] scratch_data Scratch data containing finite element
   * information at quadrature points for efficient assembly
   * @param[out] copy_data Copy data structure for storing local assembly
   * results before copying to global structures
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief Set up the vector of assembler functions for VANS equations.
   *
   * Initializes the assembler objects responsible for computing different
   * terms in the VANS equations, including particle-fluid interaction
   * assemblers and void fraction-dependent terms.
   */
  void
  setup_assemblers() override;

  /**
   * @brief Copy local matrix contributions to the global system matrix.
   *
   * Transfers the locally assembled matrix contributions from a single
   * cell to the appropriate locations in the global system matrix,
   * handling constraint applications and matrix sparsity patterns.
   *
   * @param[in] copy_data Local assembly data containing matrix contributions
   * and associated degree of freedom indices
   */
  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copy local right-hand side contributions to the global vector.
   *
   * Transfers the locally assembled right-hand side contributions from
   * a single cell to the appropriate locations in the global right-hand
   * side vector, handling constraint applications.
   *
   * @param[in] copy_data Local assembly data containing right-hand side
   * contributions and associated degree of freedom indices
   */
  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Gather additional output fields for visualization.
   *
   * Adds void fraction and particle velocity fields to the standard
   * fluid dynamics output for comprehensive visualization of the
   * coupled CFD-DEM system.
   *
   * @return Vector of OutputStruct objects containing additional fields
   * for VTU output file generation
   */
  virtual std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;


  /**
   * @brief Member Variables
   */

  /// CFD-DEM simulation parameters containing all configuration options
  CFDDEMSimulationParameters<dim> cfd_dem_simulation_parameters;

  /// Mapping object for particle position calculations and projections
  MappingQGeneric<dim> particle_mapping;

  /// Vector of assembler objects for particle-fluid interaction terms
  std::vector<std::shared_ptr<ParticleFluidAssemblerBase<dim>>>
    particle_fluid_assemblers;

  /// Flag enabling Pressure Stabilizing Petrov-Galerkin (PSPG) stabilization
  const bool PSPG = true;

  /// Flag enabling Streamline Upwind Petrov-Galerkin (SUPG) stabilization
  const bool SUPG = true;

  /// Scaling factor for GLS velocity stabilization terms
  const double GLS_u_scale = 1;

  /// Pressure drop across the computational domain
  double pressure_drop;

  /// Mapping from vertex indices to sets of cells containing each vertex
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;

  /// Mapping from vertex indices to sets of periodic cells containing each
  /// vertex
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;

  /// Particle handler for managing particle data and operations
  Particles::ParticleHandler<dim, dim> particle_handler;

  /// Particle projector for calculating void fraction and particle effects
  ParticleProjector<dim> particle_projector;

  /// Flag indicating whether the domain has periodic boundary conditions
  bool has_periodic_boundaries;

  /// Offset vector for periodic boundary condition calculations
  Tensor<1, dim> periodic_offset;

  /// Direction index for periodic boundary conditions
  unsigned int periodic_direction;
};

#endif
