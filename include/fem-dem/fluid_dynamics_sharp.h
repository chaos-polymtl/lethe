// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_gls_sharp_h
#define lethe_gls_sharp_h

#include <core/ib_particle.h>
#include <core/ib_stencil.h>
#include <core/shape.h>
#include <core/vector.h>

#include <solvers/fluid_dynamics_matrix_based.h>

#include <dem/data_containers.h>
#include <dem/log_collision_data.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/visualization.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/ib_particles_dem.h>
#include <fem-dem/sharp_ib_periodic_boundaries.h>

#include <deal.II/particles/particle_handler.h>

#include <unordered_set>

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization and
 * Sharp-Edge immersed boundaries
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 */

template <int dim>
class FluidDynamicsSharp : public FluidDynamicsMatrixBased<dim>
{
public:
  FluidDynamicsSharp(CFDDEMSimulationParameters<dim> &nsparam);

  ~FluidDynamicsSharp();

  void
  solve() override;

  /**
   * @brief Call for the assembly of the matrix
   */
  void
  assemble_system_matrix() override
  {
    assemble_matrix_and_rhs();
  }

  /**
   * @brief Call for the assembly of the right-hand side
   */
  void
  assemble_system_rhs() override
  {
    assemble_rhs();
    sharp_edge();
  }

protected:
  /**
   * @brief Override to extend the sparsity pattern for periodic IB coupling.
   *
   * The base class builds the sparsity pattern with
   * keep_constrained_dofs=false, which omits constrained (periodic slave) DoF
   * rows entirely. sharp_edge() re-imposes constraint equations on slave rows
   * for cut cells near the periodic boundary, writing two types of entry:
   *
   *   A(slave, slave) += sum_line           — diagonal identity scaling
   *   A(slave, master) += weight * sum_line  — off-diagonal coupling term
   *
   * Both must be in the sparsity pattern. We add them here for every
   * constrained DoF, since cut-cell membership is not yet known at sparsity
   * setup time. We do NOT use keep_constrained_dofs=true because that adds
   * the full element coupling for all hanging-node slave rows, which would
   * massively inflate the sparsity pattern.
   */
  void
  setup_dofs_fd() override;

private:
  struct ParticlePeriodicFrame;
  struct CutCellInfo;
  struct OverconstrainedFluidCellInfo;
  struct InsideCellInfo;
  struct PeriodicFrameDofCacheKey;

  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief sets up the vector of assembler functions
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle. Refer to 2 different assembler
   * depending on what type of equations is assembled inside the particles.
   */
  void
  setup_assemblers() override;


  /**
   * @brief Copy local cell information to global matrix
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copy local cell rhs information to global rhs
   * This function is modified compared to the GLS function to take into account
   * the cells that are cut or inside a particle
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Call for the assembly of the matrix and the right hand side
   *
   * @deprecated This function is to be deprecated when the non-linear solvers
   * have been refactored to call for rhs and matrix assembly seperately.
   */

  /*
   Modified version of assemble_matrix_and_rhs to include the presence of
   extra steps. For more detail see the same function in the
   fluid_dynamics_matrix_based.h solver.
   */

  virtual void
  assemble_matrix_and_rhs()
  {
    if (some_particles_are_coupled)
      {
        if (this->get_current_newton_iteration() != 0)
          {
            force_on_ib();
          }
        integrate_particles();
        update_precalculations_for_ib();
      }
    // Do the cut cell mapping only if it is the first Newton iteration or
    // if the explicit position integration is false.
    if (this->get_current_newton_iteration() == 0 ||
        this->simulation_parameters.particlesParameters
            ->explicit_position_integration_calculation == false)
      {
        if (all_spheres)
          optimized_generate_cut_cells_map();
        else
          generate_cut_cells_map();
      }
    this->FluidDynamicsMatrixBased<dim>::assemble_system_matrix();

    sharp_edge();
  }


  /**
   * @brief Call for the assembly of the right hand side
   *
   * @deprecated This function is to be deprecated when the non-linear solvers
   * have been refactored to call for rhs and matrix assembly seperately.
   *
   * Modified version of assemble_matrix_and_rhs to include the presence of
   * extra steps. For more detail see the same function in the
   * fluid_dynamics_matrix_based.h solver.
   */
  virtual void
  assemble_rhs()
  {
    this->FluidDynamicsMatrixBased<dim>::assemble_system_rhs();
  }

  /**
   * @brief Build the vertex-to-cell map used by Sharp-IB cell searches.
   */
  void
  vertices_cell_mapping();

  /**
   * @brief Instantiate Sharp-IB particles from the parameter file.
   */
  void
  define_particles();

  /**
   * @brief Evaluate hydrodynamic forces and torques on immersed-boundary
   * particles.
   */
  void
  force_on_ib();


  /**
   * @brief Impose Sharp-IB constraints on the assembled fluid system.
   *
   * The details of this approach are presented in: L. Barbeau, S.
   * Étienne, C. Béguin & B. Blais, «Development of a high-order continuous
   * Galerkin sharp-interface immersed boundary method and its application to
   * incompressible flow problems,» Computers & Fluids, 2020, in press, ref.
   * CAF-D-20-00773
   */
  void
  sharp_edge();

  /**
   * @brief Write per-particle force and motion history for Sharp-IB output.
   */
  void
  write_force_ib();

  /**
   * @brief Advance particle translation and rotation using the current force
   * and torque state.
   */
  void
  integrate_particles();

  /**
   * @brief Finalize the current particle time step and update particle history
   * states.
   */
  void
  finish_time_step_particles();

  /**
   * @brief Mirror DEM-style particle output and post-processing for
   * Sharp-IB particles.
   */
  void
  handle_dem_particle_output_and_postprocessing();

  /**
   * @brief Build a temporary DEM particle handler from the current Sharp-IB
   * particles.
   *
   * @param particle_handler Particle handler that stores the locally owned and
   * ghosted Sharp-IB particles.
   * @param particle_container Map from particle ids to particle iterators in
   * the temporary particle handler.
   * @param local_particle_ids Set of locally owned particle ids on the current
   * MPI rank.
   */
  void
  build_ib_particle_handler(
    Particles::ParticleHandler<dim, dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                     &particle_container,
    std::unordered_set<unsigned int> &local_particle_ids) const;

  /**
   * @brief Build DEM-style contact containers for Sharp-IB particles.
   *
   * @param particle_container Map from particle ids to particle iterators in
   * the temporary particle handler.
   * @param local_particle_ids Set of locally owned particle ids on the current
   * MPI rank.
   * @param local_adjacent_particles Particle-particle contacts where both
   * particles are locally owned.
   * @param ghost_adjacent_particles Particle-particle contacts involving a
   * ghost particle.
   * @param particle_wall_in_contact Particle-wall contacts for locally owned
   * particles.
   */
  void
  build_contact_containers(
    const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                           &particle_container,
    const std::unordered_set<unsigned int> &local_particle_ids,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
      &particle_wall_in_contact) const;

  /**
   * @brief Fill DEM particle properties from a Sharp-IB particle state.
   *
   * @param particle Sharp-IB particle whose properties are exported.
   * @param properties DEM property vector to populate.
   */
  void
  fill_particle_properties(const IBParticle<dim> &particle,
                           std::vector<double>   &properties) const;

  /**
   * @brief Evaluate the L2 error on the Sharp-IB computational domain.
   *
   * This function differs from the standard GLS counterpart because cells cut
   * by immersed boundaries are excluded from the error evaluation. See the
   * "computational domain" definition in: L. Barbeau, S. Étienne, C. Béguin &
   * B. Blais, «Development of a high-order continuous Galerkin sharp-interface
   * immersed boundary method and its application to incompressible flow
   * problems,» Computers & Fluids, 2020, in press, ref. CAF-D-20-0077
   *
   * @return Pair of L2 errors for velocity and pressure.
   */
  std::pair<double, double>
  calculate_L2_error_particles();


  /**
   * @brief Post-process Sharp-IB results using the Sharp-specific error and
   * output paths.
   */
  virtual void
  postprocess_fd(bool firstIter) override;


  /**
   * @brief Adds levelset to output files.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  virtual std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;

  /**
   * @brief Adapt the mesh around each immersed-boundary particle.
   *
   * The refinement zone is defined by a ring in 2D and a shell in 3D. The
   * outside and inside radius of this ring\shell are defined relative to the
   * particle radius by the immersed-boundary parameters
   * "refine mesh inside radius factor" and
   * "refine mesh outside radius factor". When enabled, distance-based
   * coarsening is also applied outside the corresponding coarsening threshold.
   *
   * @param initial_refinement Whether this is the initial
   * refinement cycle.
   */
  void
  mesh_adapt_ib(const bool initial_refinement);


  /**
   * @brief Check whether every Sharp-IB particle is spherical.
   *
   * If all particles are spheres, the optimized cut-cell path can be used.
   * Otherwise, Sharp falls back to the generic cut-cell classification.
   */
  void
  check_whether_all_particles_are_sphere();

  /**
   * @brief Build cut-cell, inside-cell, and overconstrained-cell maps for the
   * current particle state.
   */
  void
  generate_cut_cells_map();

  /**
   * @brief Build the optimized cut-cell maps used when all particles are
   * spheres.
   */
  void
  optimized_generate_cut_cells_map();

  /**
   * @brief Determine whether a cell is a candidate inside-cell or cut-cell for
   * a particle.
   *
   * @param cell Cell under consideration.
   * @param p_id Particle index.
   * @return Pair of flags `(cell_is_inside, cell_is_cut)`.
   */
  std::pair<bool, bool>
  generate_cut_cell_candidates(
    const typename DoFHandler<dim>::cell_iterator &cell,
    const unsigned int                             p_id);

  /**
   * @brief Check whether a cell is cut by an immersed-boundary particle.
   *
   * @param cell Cell under consideration.
   * @param local_dof_indices Local DoF indices for @p cell.
   * @param support_points Support-point map for the flow DoFs.
   * @return Tuple `(cell_is_cut, particle_id, local_dof_indices)`.
   */
  std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
  cell_cut(const typename DoFHandler<dim>::active_cell_iterator &cell,
           std::vector<types::global_dof_index>          &local_dof_indices,
           std::map<types::global_dof_index, Point<dim>> &support_points);

  /**
   * @brief Check whether a cell is cut by a particle whose geometry does not
   * provide a signed level-set test.
   *
   * @param cell Cell under consideration.
   * @param support_points Support-point map for the flow DoFs.
   * @param p Particle index used in the check.
   * @return `true` if the cell is cut by particle @p p.
   */
  bool
  cell_cut_by_p_absolute_distance(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    std::map<types::global_dof_index, Point<dim>>        &support_points,
    unsigned int                                          p);

  /**
   * @brief Check whether a cell lies fully inside an immersed-boundary
   * particle.
   *
   * @param cell Cell under consideration.
   * @param local_dof_indices Local DoF indices for @p cell.
   * @param support_points Support-point map for the flow DoFs.
   * @return Tuple `(cell_is_inside, particle_id, local_dof_indices)`.
   */
  std::tuple<bool, unsigned int, std::vector<types::global_dof_index>>
  cell_inside(const typename DoFHandler<dim>::active_cell_iterator &cell,
              std::vector<types::global_dof_index>          &local_dof_indices,
              std::map<types::global_dof_index, Point<dim>> &support_points);



  /**
   * @brief Write a Sharp-IB checkpoint that supports simulation restart.
   *
   * This function stores all previous particle states in one file. Each
   * row corresponds to one particle state. The file is structured as follows:
   *
   *
   * P0 state at time t
   *
   * P0 state at time t-dt
   *
   * P0 state at time t-2dt
   *
   * P1 state at time t
   *
   * P1 state at time t-dt
   *
   * P1 state at time t-2dt
   *
   * etc
   */
  virtual void
  write_checkpoint() override;

  /**
   * @brief Read a Sharp-IB checkpoint and initialize a restart.
   *
   * See write_checkpoint() for the file structure.
   */
  virtual void
  read_checkpoint() override;


  /**
   * @brief Read a particle file and populate the initial Sharp-IB particle
   * state.
   *
   * The file must contain the following information for each particle (the
   * header must be defined accordingly):
   * type shape_argument_0 shape_argument_1 shape_argument_2 p_x p_y p_z v_x v_y
   * v_z omega_x omega_y omega_z orientation_x orientation_y orientation_z
   * density inertia pressure_x pressure_y pressure_z youngs_modulus
   * restitution_coefficient friction_coefficient poisson_ratio
   * rolling_friction_coefficient
   */
  void
  load_particles_from_file();

  /**
   * @brief Regroup and organize the refinement process around the IB particle.
   *
   * @param initial_refinement Whether this call performs the initial
   * refinement before time stepping starts.
   */
  void
  refinement_control(const bool initial_refinement);


  /**
   * @brief Check whether a point lies inside a given active cell.
   *
   * @param cell Active cell under consideration.
   * @param point Point to test.
   * @return `true` if @p point lies inside @p cell.
   */
  bool
  point_inside_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                    Point<dim> point);

  /**
   * @brief Override the nonlinear residual to include particle-coupling
   * residuals.
   *
   * @return The larger of the fluid residual and the particle residual scaled
   * to the fluid nonlinear tolerance.
   */
  double
  get_current_residual() override
  {
    double scaling = this->simulation_parameters.physics_solving_strategy
                       .at(PhysicsID::fluid_dynamics)
                       .tolerance /
                     this->simulation_parameters.particlesParameters
                       ->particle_nonlinear_tolerance;
    return std::max(this->system_rhs.l2_norm(), particle_residual * scaling);
  }

  /**
   * @brief Update cached geometric data for every immersed particle.
   */
  void
  update_precalculations_for_ib();

  /**
   * @brief Validate and cache the periodic configuration used by Sharp-IB.
   *
   * This phase depends only on DEM/CFD boundary-condition definitions. It must
   * run before setup_dofs_fd() so the sparsity pattern knows whether periodic
   * Sharp-specific entries are required.
   */
  void
  initialize_sharp_ib_periodic_configuration();

  /**
   * @brief Update periodic domain bounds and offset for the active mesh.
   *
   * This phase depends on the current triangulation and must be called after
   * the mesh exists in its current refinement state.
   */
  void
  update_sharp_ib_periodic_geometry();

  /**
   * @brief Add matrix entry, expanding column constraints for periodic DoFs.
   *
   * When col_dof is a periodic-constrained (slave) DoF, the sparsity pattern
   * (built with keep_constrained_dofs=false) does not include (row, slave)
   * entries. This function expands the column constraint to redirect the
   * entry to (row, master_col), which IS in the sparsity pattern.
   *
   * Row constraints are NOT expanded. sharp_edge() deliberately writes
   * equations on constrained DoF rows; redirecting to master rows would
   * corrupt master DoF equations. After the solve, constraints.distribute()
   * correctly sets constrained DoF values from their masters.
   *
   * @param row_dof Global DoF index for matrix row
   * @param col_dof Global DoF index for matrix column
   * @param value Value to add to matrix entry
   */
  void
  add_matrix_entry_with_periodic_expansion(
    const types::global_dof_index row_dof,
    const types::global_dof_index col_dof,
    const double                  value);

  /**
   * @brief Get periodic images of a particle for refinement/coupling.
   *
   * Given a particle position and radius, returns a list of periodic image
   * positions that may affect cells across the periodic boundary.
   *
   * @param particle_position Center position of particle
   * @param particle_radius Radius of particle (or support radius for effect)
   * @return Vector of periodic image positions (empty if no images needed)
   */
  std::vector<Point<dim>>
  get_periodic_particle_images(const Point<dim> &particle_position,
                               const double      particle_radius) const;

  /**
   * @brief Get the periodic shifts needed to evaluate all relevant particle
   * frames for a given support region.
   *
   * The returned list always begins with the zero shift corresponding to the
   * physical (primary) particle position. Callers must not prepend or skip
   * this zero shift: iterating the full returned list covers both the primary
   * frame and any periodic images that overlap the support region.
   *
   * When periodic boundaries are disabled or the support region does not cross
   * the periodic boundary, the list contains only the zero shift.
   *
   * @param particle_position Particle centroid in the primary domain.
   * @param support_radius Radius used to detect whether a periodic image can
   * influence the current operation.
   * @return List of periodic shifts, starting with the zero shift for the
   * primary frame.
   */
  std::vector<Tensor<1, dim>>
  get_periodic_particle_shifts(const Point<dim> &particle_position,
                               const double      support_radius) const;

  /**
   * @brief Convert a periodic particle shift into a discrete cache-frame
   * index.
   *
   * Sharp-IB currently supports a single periodic pair, so every particle
   * evaluation frame is either the primary frame (`0`) or one wrapped image on
   * the positive/negative side of the periodic direction (`+1` or `-1`).
   * Cache keys use this discrete index instead of the raw shift tensor so
   * stencil/force reuse remains stable and unambiguous.
   *
   * @param periodic_shift Particle-frame shift relative to the primary
   * position.
   * @return `0` for the primary frame, `+1` or `-1` for wrapped periodic
   * images.
   */
  int
  get_periodic_frame_index(const Tensor<1, dim> &periodic_shift) const;

  /**
   * @brief Apply a periodic shift to a particle position for temporary
   * geometric evaluations.
   *
   * This function mutates the particle's position. Every call must be paired
   * with a matching restore of @p primary_particle_position before the
   * particle is used for any purpose other than the immediate evaluation.
   * Failing to restore leaves the particle in a shifted frame for all
   * subsequent operations in the current time step.
   *
   * @param particle_id Particle index.
   * @param primary_particle_position Particle centroid in the primary domain.
   * @param periodic_shift Shift relative to @p primary_particle_position.
   * Pass a zero-initialized tensor to restore the particle to its primary
   * frame position.
   */
  void
  set_particle_position_from_periodic_shift(
    const unsigned int    particle_id,
    const Point<dim>     &primary_particle_position,
    const Tensor<1, dim> &periodic_shift);

  /**
   * @brief Register the periodic frame that cut a cell.
   *
   * Sharp-IB currently supports exactly one periodic pair. A cell should
   * therefore be cut by at most one periodic frame of a given particle.
   * Storing the frame explicitly here keeps later assembly/force paths from
   * re-guessing which wrapped particle image generated the cut-cell data.
   *
   * @param cut_cell_info Cut-cell metadata to update.
   * @param particle_id Particle index.
   * @param periodic_shift Shift of the periodic frame that cut the cell.
   */
  void
  register_cutting_particle_frame(CutCellInfo          &cut_cell_info,
                                  const unsigned int    particle_id,
                                  const Tensor<1, dim> &periodic_shift);

  /**
   * @brief Recover the stored periodic frame for a particle that cuts a cell.
   *
   * @param cut_cell_info Cut-cell metadata for the cell under consideration.
   * @param particle_id Particle index.
   * @return Periodic shift used when that particle was classified as cutting
   * the cell.
   */
  Tensor<1, dim>
  get_periodic_shift_for_particle_frame(const CutCellInfo &cut_cell_info,
                                        const unsigned int particle_id) const;

  /**
   * @brief Defines a struct with methods that allow the generation of a visualisation of the IB_particles. This is equivalent to the corresponding class in the DEM solver.
   */
  struct Visualization_IB : public dealii::DataOutInterface<0, dim>
  {
  public:
    /**
     * Carries out building the patches of properties of particles for
     * visualization
     *
     * @param particles The vector fo IB_particles
     */
    void
    build_patches(std::vector<IBParticle<dim>> particles);


    ~Visualization_IB();

  private:
    /**
     * Implementation of the corresponding function of the base class.
     */
    virtual const std::vector<DataOutBase::Patch<0, dim>> &
    get_patches() const override;

    /**
     * Implementation of the corresponding function of the base class.
     */
    virtual std::vector<std::string>
    get_dataset_names() const override;

    virtual std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
    get_nonscalar_data_ranges() const override;



    /**
     * Output information that is filled by build_patches() and
     * written by the write function of the base class.
     */
    std::vector<DataOutBase::Patch<0, dim>> patches;

    /**
     * A list of field names for all data components stored in patches.
     */
    std::vector<std::string> dataset_names;

    /**
     * Store which of the data fields are vectors.
     */

    std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      vector_datasets;


    /**
     * Particle properties that are written in output files
     */
    std::vector<std::pair<std::string, int>> properties_to_write;
  };



  /**
   * Members
   */

private:
  struct ParticlePeriodicFrame
  {
    unsigned int   particle_id = 0;
    Tensor<1, dim> periodic_shift;
  };

  struct CutCellInfo
  {
    bool                               is_cut              = false;
    unsigned int                       primary_particle_id = 0;
    std::vector<ParticlePeriodicFrame> cutting_particle_frames;
  };

  struct OverconstrainedFluidCellInfo
  {
    bool           is_overconstrained  = false;
    unsigned int   particle_id         = 0;
    double         distance_to_surface = 0.0;
    Tensor<1, dim> periodic_shift;
  };

  struct InsideCellInfo
  {
    bool           is_inside   = false;
    unsigned int   particle_id = 0;
    Tensor<1, dim> periodic_shift;
  };

  struct PeriodicFrameDofCacheKey
  {
    types::global_dof_index dof_index   = 0;
    unsigned int            particle_id = numbers::invalid_unsigned_int;
    int                     periodic_frame_index = 0;

    bool
    operator<(const PeriodicFrameDofCacheKey &other) const
    {
      if (dof_index != other.dof_index)
        return dof_index < other.dof_index;
      if (particle_id != other.particle_id)
        return particle_id < other.particle_id;
      return periodic_frame_index < other.periodic_frame_index;
    }
  };

  // Parameters
  CFDDEMSimulationParameters<dim> cfd_dem_parameters;

  bool all_spheres;

  // Bool that check if some particle are coupled.
  bool some_particles_are_coupled;

  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;

  // For each active cell, store whether it is cut by any immersed particle,
  // which particle owns the primary cut-cell treatment, and which periodic
  // particle frames generated the cut classification.
  std::map<typename DoFHandler<dim>::active_cell_iterator, CutCellInfo>
    cut_cells_map;

  // For each active cell, store whether the fluid cell is overconstrained, the
  // particle responsible for that state, the distance to the nearest particle
  // surface, and the periodic frame in which that classification was made.
  std::map<typename DoFHandler<dim>::active_cell_iterator,
           OverconstrainedFluidCellInfo>
    overconstrained_fluid_cell_map;

  /*
   * These vectors and map are used to keep track of the DOFs that are
   * overconstrained
   */
  GlobalVectorType             local_dof_overconstrained;
  GlobalVectorType             dof_overconstrained;
  std::map<unsigned int, bool> dof_with_more_then_one_particle;

  // For each active cell, store whether it lies fully inside an immersed
  // particle, which particle owns it, and the periodic frame used for the
  // inside-cell classification.
  std::map<typename DoFHandler<dim>::active_cell_iterator, InsideCellInfo>
    cells_inside_map;
  /*
   * Cache of stencil cells already located for a specific
   * (DoF, particle, periodic frame) tuple during sharp_edge().
   *
   * The periodic frame is part of the key because the same DoF can be reached
   * from different wrapped images of the same particle near the periodic
   * interface, and those frames can require different stencil cells.
   */
  std::map<PeriodicFrameDofCacheKey,
           typename DoFHandler<dim>::active_cell_iterator>
    ib_done;

  // Special assembler of the cells inside an IB particle
  std::vector<std::shared_ptr<NavierStokesAssemblerBase<dim>>>
    assemblers_inside_ib;

  PVDHandler ib_particles_pvdhandler;
  // DEM-style post-processing handlers rebuilt from Sharp-IB particles.
  PVDHandler ib_particles_pvdhandler_force_chains;

  std::vector<IBParticle<dim>> particles;
  double                       particle_residual;

  std::vector<TableHandler> table_p;
  TableHandler              table_all_p;

  // Object used to sub-time step the particle dynamics to allow contact between
  // particles.
  IBParticlesDEM<dim> ib_dem;

  // DEM-style collision data used by Sharp post-processing.
  OngoingCollisionLog<dim>   ongoing_collision_log;
  CompletedCollisionLog<dim> collision_event_log;

  // Function that describes all solids signed distance functions together
  std::shared_ptr<Shape<dim>> combined_shapes;

  // Postprocessors to output the signed distance function of the immersed
  // solids
  std::shared_ptr<LevelsetPostprocessor<dim>> levelset_postprocessor;
  std::shared_ptr<LevelsetGradientPostprocessor<dim>>
    levelset_gradient_postprocessor;

  // Combined periodic boundary handler for sharp IB coupling
  SharpIBPeriodicBoundaries<dim> sharp_ib_periodic_boundaries;

  // Cache for domain bounds in periodic direction (used for wrapping)
  double periodic_domain_lower;
  double periodic_domain_upper;
};


#endif
