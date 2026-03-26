// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_h
#define lethe_cls_h

#include <core/bdf.h>
#include <core/interface_tools.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/cls_assemblers.h>
#include <solvers/cls_filter.h>
#include <solvers/cls_linear_subequations_solver.h>
#include <solvers/cls_scratch_data.h>
#include <solvers/cls_subequations_interface.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/signed_distance_transformation.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>


DeclException1(
  InvalidNumberOfFluid,
  int,
  << "The CLS physics is enabled, but the number of fluids is set to " << arg1
  << ". The CLS solver only supports 2 fluids.");

DeclException1(
  CLSBoundaryConditionMissing,
  types::boundary_id,
  << "The boundary id: " << arg1
  << " is defined in the triangulation, but not as a boundary condition for the CLS physics. Lethe does not assign a default boundary condition to boundary ids. Every boundary id defined within the triangulation must have a corresponding boundary condition defined in the input file.");

DeclExceptionMsg(
  UnsupportedRegularization,
  "The CLS physics has been set to use DG and the latter implementation currently does not support any interface reinitialization mechanism.");

DeclExceptionMsg(
  UnsupportedRegularizationWithSimplex,
  "The CLS physics has been set to use simplex and the latter implementation currently does not support the geometric interface reinitialization mechanism.");

DeclExceptionMsg(
  UnsupportedInitialProjection,
  "The CLS physics has been set to use DG and the latter implementation currently does not support defining an initial condition with a projection.");

template <int dim>
class ConservativeLevelSet : public AuxiliaryPhysics<dim, GlobalVectorType>
{
public:
  /**
   * @brief CLS - Base constructor.
   */
  ConservativeLevelSet(
    MultiphysicsInterface<dim>      *multiphysics_interface,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control);

  /**
   * @brief CLS - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  ~ConservativeLevelSet() = default;

  /**
   * @brief Call for the assembly of the right-hand side
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver
   */
  void
  assemble_rhs();

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  virtual std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;

  /**
   * @brief Calculates the L2 error of the solution
   */
  double
  calculate_L2_error();

  /**
   * @brief Calculates the volume and mass for a given fluid phase.
   * Used for conservation monitoring.
   *
   * @param solution CLS solution (phase indicator)
   *
   * @param current_solution_fd current solution for the fluid dynamics
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1) corresponding to
   * the phase of interest.
   */
  template <typename VectorType>
  void
  calculate_volume_and_mass(const GlobalVectorType &solution,
                            const VectorType       &current_solution_fd,
                            const Parameters::FluidIndicator monitored_fluid);

  /**
   * @brief Calculates the momentum for a given fluid phase.
   * Used for conservation monitoring.
   *
   * @param[in] solution CLS solution (phase indicator)
   *
   * @param[in] current_solution_fd current solution for the fluid dynamics
   *
   * @param[in] monitored_fluid Fluid indicator (fluid0 or fluid1) corresponding
   * to the phase of interest.
   *
   * @return A tensor<1,dim> corresponding to the entry_string in the prm file.
   */
  template <typename VectorType>
  Tensor<1, dim>
  calculate_momentum(const GlobalVectorType          &solution,
                     const VectorType                &current_solution_fd,
                     const Parameters::FluidIndicator monitored_fluid);

  /**
   * @brief Calculates the barycenter of the fluid and its velocity
   *
   * @param[in] solution CLS solution
   *
   * @param[in] current_solution_fd Fluid dynamics solution
   *
   */
  template <typename VectorType>
  std::pair<Tensor<1, dim>, Tensor<1, dim>>
  calculate_barycenter(const GlobalVectorType &solution,
                       const VectorType       &current_solution_fd);

  /**
   * @brief Computes the capillary time-step constraint
   *
   * The time-step constraint for capillary-driven multiphase flows
   * is computed according to [1]:
   * \f[
   * \Delta t_\sigma = \sqrt{\frac{(\rho_0 + \rho_1) h^3}{2 \pi \sigma}}
   * \f]
   *
   * @param[in] pressure_dependent Boolean indicating if the pressure field is
   * required. The DensityIsothermalIdealGas model requires the pressure field
   * to compute the density.
   *
   * @param[in] temperature_dependent Boolean indicating if the temperature
   * field is required. The SurfaceTensionLinear and SurfaceTensionPhaseChange
   * models depend on the temperature field.
   *
   * @return Computed value of the capillary time-step constraint
   *
   * Reference:
   *
   * [1] F. Denner, F. Evrard, and B. van Wachem, “Breaching the capillary
   * time-step constraint using a coupled VOF method with implicit surface
   * tension,” J. Comput. Phys., vol. 459, p. 111128, Jun. 2022,
   * doi: 10.1016/j.jcp.2022.111128
   */
  double
  compute_capillary_time_step_constraint(const bool pressure_dependent,
                                         const bool temperature_dependent);

  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  void
  finish_simulation() override;

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Carry out modifications on the auxiliary physic solution.
   */
  void
  modify_solution() override;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  void
  postprocess(bool first_iteration) override;


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  void
  pre_mesh_adaptation() override;

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  void
  post_mesh_adaptation() override;


  /**
   * @brief Compute the error estimator for mesh refinement.
   *
   * @param[in] ivar The current element of the map
   * simulation_parameters.mesh_adaptation.variables
   *
   * @param[in,out] estimated_error_per_cell The deal.II vector of
   * estimated_error_per_cell
   */
  void
  compute_error_estimate(
    const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                          &ivar,
    dealii::Vector<float> &estimated_error_per_cell) override;

  /**
   * @brief Compute the Kelly error estimator on the phase variable for mesh refinement.
   * See :
   * https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html
   * for more information on the Kelly error estimator.
   *
   * @param[in,out] estimated_error_per_cell The deal.II vector of
   * estimated_error_per_cell
   * @param[in] component_mask The component mask corresponding to the phase
   * variable
   */
  void
  compute_kelly(dealii::Vector<float> &estimated_error_per_cell,
                const ComponentMask   &component_mask);

  /**
   * @brief Prepares auxiliary physics to write checkpoint
   */
  void
  write_checkpoint() override;


  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  void
  read_checkpoint() override;

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for the CLS solver.
   *
   * @return Structure containing a vector of references to TableHandler objects
   * that needs to be serialized/deserialized for the CLS solver, and their
   * corresponding file names.
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  void
  setup_dofs() override;

  /**
   * @brief  defined the zero constraints used to solve the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief  defined the non zero constraints used to solve the problem.
   */
  void
  define_non_zero_constraints();

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionally support
   * the use of L2 projection or steady-state solutions.
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Update non zero constraints if the boundary is time dependent
   */
  void
  update_boundary_conditions() override;

  /**
   * @brief Call for the solution of the linear system of equation using a strategy appropriate
   * to the auxiliary physics
   */
  void
  solve_linear_system() override;

  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   * NB : dof_handler and present_solution are passed to the multiphysics
   * interface at the end of the setup_dofs method
   */

  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return *dof_handler;
  }

  GlobalVectorType &
  get_evaluation_point() override
  {
    return evaluation_point;
  }

  GlobalVectorType &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  }

  GlobalVectorType &
  get_newton_update() override
  {
    return newton_update;
  }

  GlobalVectorType &
  get_present_solution() override
  {
    return *present_solution;
  }

  GlobalVectorType &
  get_system_rhs() override
  {
    return system_rhs;
  }

  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }

  const DoFHandler<dim> &
  get_projected_phase_indicator_gradient_dof_handler()
  {
    return this->cls_subequations_interface->get_dof_handler(
      CLSSubequationsID::phase_gradient_projection);
  }

  const DoFHandler<dim> &
  get_curvature_dof_handler()
  {
    return this->cls_subequations_interface->get_dof_handler(
      CLSSubequationsID::curvature_projection);
  }

  const GlobalVectorType &
  get_projected_phase_indicator_gradient_solution()
  {
    return this->cls_subequations_interface->get_solution(
      CLSSubequationsID::phase_gradient_projection);
  }

  const GlobalVectorType &
  get_curvature_solution()
  {
    return this->cls_subequations_interface->get_solution(
      CLSSubequationsID::curvature_projection);
  }
  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   */
  void
  output_newton_update_norms() override
  {
    this->pcout << "\t  ||dphi||_L2 = " << std::setw(6)
                << newton_update.l2_norm() << std::setw(6)
                << "\t  ||dphi||_Linfty = " << newton_update.linfty_norm()
                << std::endl;
  }

  /**
   * @brief Return the metric for residual rescaling. By default, should return 1.
   * If the rescale_residual_by_volume is set to true, the method
   * returns the square root of the global volume of the triangulation.
   *
   * @return Rescale metric.
   */
  double
  get_residual_rescale_metric() const override
  {
    return simulation_parameters.linear_solver.at(PhysicsID::CLS)
               .rescale_residual_by_volume ?
             std::sqrt(
               GridTools::volume(*this->triangulation, *this->mapping)) :
             1.;
  }

private:
  /**
   * @brief Verify consistency of the input parameters for boundary
   * conditions to ensure that for every boundary condition within the
   * triangulation, a boundary condition has been specified in the input file.
   */
  void
  verify_consistency_of_boundary_conditions()
  {
    // Sanity check all of the boundary conditions of the triangulation to
    // ensure that they have a type.
    std::vector<types::boundary_id> boundary_ids_in_triangulation =
      this->triangulation->get_boundary_ids();
    for (auto const &boundary_id_in_tria : boundary_ids_in_triangulation)
      {
        AssertThrow(simulation_parameters.boundary_conditions_cls.type.find(
                      boundary_id_in_tria) !=
                      simulation_parameters.boundary_conditions_cls.type.end(),
                    CLSBoundaryConditionMissing(boundary_id_in_tria));
      }
  }

  /**
   *  @brief Assembles the matrix associated with the solver.
   */
  void
  assemble_system_matrix() override;

  /**
   *  @brief Assemble the matrix associated with the solver when CG elements are used.
   */
  void
  assemble_system_matrix_cg();

  /**
   *  @brief Assemble the matrix associated with the solver when DG elements are used.
   */
  void
  assemble_system_matrix_dg();

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   *  @brief Assemble the rhs associated with the solver when CG elements are used.
   */
  void
  assemble_system_rhs_cg();

  /**
   *  @brief Assemble the rhs associated with the solver when DG elements are used.
   */
  void
  assemble_system_rhs_dg();

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
   * See the documentation for CLSScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    CLSScratchData<dim>                                  &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for CLSScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    CLSScratchData<dim>                                  &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief sets up the vector of assembler functions
   */
  virtual void
  setup_assemblers();


  /**
   * @brief Copy local cell information to global matrix
   */

  virtual void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Copy local cell rhs information to global rhs
   */

  virtual void
  copy_local_rhs_to_global_rhs(const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Limit the phase indicators between 0 and 1. This is necessary before interface sharpening.
   * More information can be found in step_41 of deal.II tutorials:
   * https://www.dealii.org/current/doxygen/deal.II/step_41.html
   */
  void
  update_solution_and_constraints(GlobalVectorType &solution);

  /**
   * @brief Assemble the system for interface sharpening
   *  * This function assembles the weak form of:
   * \f$ \Phi = c ^ {(1 - \alpha)} * (\phi ^ \alpha)\f$  if \f$ 0 <=
   * \phi <= c  \ \f$
   * \f$ \Phi = 1 - (1 - c) ^ {(1 - \alpha)} * (1 - \phi) ^ \alpha \f$ if \f$c <
   * \phi <= 1 \f$
   * Reference for sharpening method
   * https://www.sciencedirect.com/science/article/pii/S0045782500002000
   *
   * @param solution CLS solution (phase indicator)
   *
   * @param sharpening_threshold Interface sharpening threshold that represents
   * the mass conservation level
   *
   */
  void
  assemble_L2_projection_interface_sharpening(
    GlobalVectorType &solution,
    const double      sharpening_threshold);

  /**
   * @brief Solves the assembled system to sharpen the interface. The linear_solver_tolerance
   * is hardcoded = 1e-15, and an ILU Preconditioner is used. After solving the
   * system, this function overwrites the solution with the sharpened solution
   *
   * @param solution CLS solution (phase indicator)
   */
  void
  solve_interface_sharpening(GlobalVectorType &solution);

  /**
   * @brief Assembles a mass_matrix which is used in update_solution_and_constraints function
   * to limit the phase indicators of cells in the range of [0,1] before
   * sharpening the interface. More information can be found in step_41 of
   * deal.II tutorials:
   * https://www.dealii.org/current/doxygen/deal.II/step_41.html
   *
   * @param mass_matrix
   */
  void
  assemble_mass_matrix(TrilinosWrappers::SparseMatrix &mass_matrix);

  /**
   * @brief Carries out interface sharpening. It is called in the modify solution function.
   * Launches sharpen_interface with the possibility to ensure conservation and
   * handles output messages.
   */
  void
  handle_interface_sharpening();

  /**
   * @brief Find the sharpening threshold to ensure mass conservation of the fluid
   * monitored, as given in the prm (CLS, subsection monitoring), by binary
   * search.
   */
  double
  find_sharpening_threshold();

  /**
   * @brief Calculate the mass deviation of the monitored fluid, between the current
   * iteration and the mass at first iteration (mass_first_iteration). Used to
   * test multiple sharpening threshold in the binary search algorithm
   * (adaptive sharpening).
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1) corresponding to
   * the phase of interest.
   *
   * @param sharpening_threshold Interface sharpening threshold that represents the
   * mass conservation level
   */
  double
  calculate_mass_deviation(const Parameters::FluidIndicator monitored_fluid,
                           const double sharpening_threshold);

  /**
   * @brief Carries out interface sharpening. It is called in the modify solution function.
   *
   * @param solution CLS solution (phase indicator)
   *
   * @param sharpening_threshold Interface sharpening threshold that represents the mass conservation level
   *
   * @param sharpen_previous_solutions Boolean true if sharpening is applied on the present
   * and past solutions, false if the sharpening is applied on the given
   * solution vector only. Used to determine the sharpening threshold by binary
   * search.
   */
  void
  sharpen_interface(GlobalVectorType &solution,
                    const double      sharpening_threshold,
                    const bool        sharpen_previous_solutions);

  /**
   * @brief Carries out the smoothing phase indicator with a projection step (to avoid a staircase interface).
   */
  void
  smooth_phase_indicator(GlobalVectorType &solution);

  /**
   * @brief Assembles the matrix and rhs for calculation of a smooth phase indicator using a projection.
   *
   * @param solution CLS solution (phase indicator)
   */
  void
  assemble_projection_phase_indicator(GlobalVectorType &solution);

  /**
   * @brief Solves smooth phase indicator system.
   * @param solution CLS solution (phase indicator)
   */
  void
  solve_projection_phase_indicator(GlobalVectorType &solution);

  /**
   * @brief Apply filter on phase indicator values.
   *
   * @param[in] original_solution CLS solution vector to which the filter is to
   * be applied.
   *
   * @param[out] filtered_solution Solution vector with filtered phase indicator
   * values.
   */
  void
  apply_phase_filter(const GlobalVectorType &original_solution,
                     GlobalVectorType       &filtered_solution);

  /**
   * @brief Reinitialize the interface between fluids using the algebraic
   * approach.
   */
  void
  reinitialize_interface_with_algebraic_method();

  /**
   * @brief Compute level-set field from the phase indicator field using an
   * inverse tanh-based transformation.
   *
   * @param[in] solution Phase indicator solution field
   *
   * @param[out] level_set_solution Level-set solution field
   */
  void
  compute_level_set_from_phase_indicator(const GlobalVectorType &solution,
                                        GlobalVectorType &level_set_solution);

  /**
   * @brief Compute the phase indicator field from level-set field using a
   * tanh-based transformation.
   *
   * @param[in] level_set_solution Level-set solution field
   *
   * @param[out] phase_indicator_solution Phase indicator solution field
   */
  void
  compute_phase_indicator_from_level_set(
    const GlobalVectorType &level_set_solution,
    GlobalVectorType       &phase_indicator_solution);

  /**
   * @brief Reinitialize the interface between fluids using the geometric
   * approach.
   */
  void
  reinitialize_interface_with_geometric_method();

  /**
   * @brief Helper function to reinit the face velocity with the adequate solution.
   * This prevents code duplication throughout the CLS class. The function looks
   * at the multiphysics interface to decide if the velocity is a block velocity
   * or a regular velocity. Furthermore, it also checks if a time-averaged
   * solution is required. Otherwise, the code here would be copied four times.
   *
   * @param[in] velocity_cell the iterator of the cell where the velocity is to
   * be reinitialized.
   *
   * @param[in] face_no the face index where the velocity is to be
   * reinitialized.
   *
   * @param[in,out] scratch_data the scratch data to be used for the
   * reinitialization.
   */
  inline void
  reinit_face_velocity_with_adequate_solution(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const unsigned int                                   &face_no,
    CLSScratchData<dim>                                  &scratch_data)
  {
    if (multiphysics->fluid_dynamics_is_block())
      {
        scratch_data.reinit_face_velocity(velocity_cell,
                                          face_no,
                                          multiphysics->get_block_solution(
                                            PhysicsID::fluid_dynamics),
                                          this->simulation_parameters.ale);
      }
    else
      {
        scratch_data.reinit_face_velocity(velocity_cell,
                                          face_no,
                                          multiphysics->get_solution(
                                            PhysicsID::fluid_dynamics),
                                          this->simulation_parameters.ale);
      }
  }


  GlobalVectorType nodal_phase_indicator_owned;

  MultiphysicsInterface<dim> *multiphysics;

  TimerOutput computing_timer;

  const SimulationParameters<dim> &simulation_parameters;

  // Core elements for the CLS simulation
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl>  simulation_control;
  std::shared_ptr<DoFHandler<dim>>    dof_handler;
  std::shared_ptr<FiniteElement<dim>> fe;
  ConvergenceTable                    error_table;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  // Solution storage
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  GlobalVectorType                  evaluation_point;
  GlobalVectorType                  local_evaluation_point;
  GlobalVectorType                  newton_update;
  std::shared_ptr<GlobalVectorType> present_solution;
  GlobalVectorType                  system_rhs;
  AffineConstraints<double>         nonzero_constraints;
  AffineConstraints<double>         bounding_constraints;
  AffineConstraints<double>         zero_constraints;
  TrilinosWrappers::SparseMatrix    system_matrix;
  std::shared_ptr<GlobalVectorType> filtered_solution;

  /// Level-set field obtained from the phase indicator field using a tanh-based
  /// transformation
  GlobalVectorType level_set;

  // Previous solutions vectors
  std::shared_ptr<std::vector<GlobalVectorType>> previous_solutions;

  // Solution transfer classes
  std::shared_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;
  std::vector<SolutionTransfer<dim, GlobalVectorType>>
    previous_solutions_transfer;

  // Phase indicator matrices for interface sharpening
  TrilinosWrappers::SparseMatrix system_matrix_phase_indicator;
  TrilinosWrappers::SparseMatrix complete_system_matrix_phase_indicator;
  GlobalVectorType               system_rhs_phase_indicator;
  GlobalVectorType               complete_system_rhs_phase_indicator;
  TrilinosWrappers::SparseMatrix mass_matrix_phase_indicator;

  // For projected phase indicator gradient (pfg), projected curvature, and
  // algebraic interface reinitialization
  std::shared_ptr<CLSSubequationsInterface<dim>> cls_subequations_interface;

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Lower and upper bounds of phase indicator
  const double phase_upper_bound = 1.0;
  const double phase_lower_bound = 0.0;

  // Conservation Analysis
  TableHandler table_monitoring_cls;
  double       volume_monitored;
  double       mass_monitored;
  double       mass_first_iteration;
  double       sharpening_threshold;

  // Barycenter analysis
  TableHandler table_barycenter;

  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<CLSAssemblerBase<dim>>> assemblers;

  // Face assemblers, used only for DG methods
  std::shared_ptr<CLSAssemblerSIPG<dim>> inner_face_assembler;

  // Phase indicator filter
  std::shared_ptr<ConservativeLevelSetFilterBase> filter;

  // Signed distance solver for geometric redistanciation
  std::shared_ptr<InterfaceTools::SignedDistanceSolver<dim, GlobalVectorType>>
    signed_distance_solver;

  // Signed distance transformation function to a phase indicator
  std::shared_ptr<SignedDistanceTransformationBase>
    signed_distance_transformation;
};



#endif
