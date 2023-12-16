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
 */

#ifndef lethe_VOF_h
#define lethe_VOF_h

#include <core/bdf.h>
#include <core/simulation_control.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_filter.h>
#include <solvers/vof_scratch_data.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/error_estimator.h>

DeclException1(
  InvalidNumberOfFluid,
  int,
  << "The VOF physics is enabled, but the number of fluids is set to " << arg1
  << ". The VOF solver only supports 2 fluids.");


template <int dim>
class VolumeOfFluid
  : public AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>
{
public:
  /**
   * @brief VOF - Base constructor.
   */
  VolumeOfFluid(MultiphysicsInterface<dim>      *multiphysics_interface,
                const SimulationParameters<dim> &p_simulation_parameters,
                std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                                   p_triangulation,
                std::shared_ptr<SimulationControl> p_simulation_control)
    : AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>(
        p_simulation_parameters.non_linear_solver.at(PhysicsID::VOF))
    , multiphysics(multiphysics_interface)
    , computing_timer(p_triangulation->get_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
    , projected_phase_fraction_gradient_dof_handler(*triangulation)
    , curvature_dof_handler(*triangulation)
    , sharpening_threshold(
        simulation_parameters.multiphysics.vof_parameters.sharpening.threshold)
  {
    AssertThrow(simulation_parameters.physical_properties_manager
                    .get_number_of_fluids() == 2,
                InvalidNumberOfFluid(
                  simulation_parameters.physical_properties_manager
                    .get_number_of_fluids()));


    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe = std::make_shared<FE_SimplexP<dim>>(
          simulation_parameters.fem_parameters.VOF_order);
        mapping         = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
        face_quadrature =
          std::make_shared<QGaussSimplex<dim - 1>>(fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe = std::make_shared<FE_Q<dim>>(
          simulation_parameters.fem_parameters.VOF_order);
        mapping = std::make_shared<MappingQ<dim>>(
          fe->degree, simulation_parameters.fem_parameters.qmapping_all);
        fe_projected_phase_fraction_gradient =
          std::make_shared<FESystem<dim>>(FE_Q<dim>(fe->degree), dim);
        fe_curvature = std::make_shared<FE_Q<dim>>(fe->degree);
        projected_phase_fraction_gradient_mapping =
          std::make_shared<MappingQ<dim>>(
            fe_projected_phase_fraction_gradient->degree,
            simulation_parameters.fem_parameters.qmapping_all);
        curvature_mapping = std::make_shared<MappingQ<dim>>(
          fe_curvature->degree,
          simulation_parameters.fem_parameters.qmapping_all);
        cell_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
        face_quadrature = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
      }

    // Allocate solution transfer
    solution_transfer =
      std::make_shared<parallel::distributed::
                         SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>(
        dof_handler);

    // Set size of previous solutions using BDF schemes information
    previous_solutions.resize(maximum_number_of_previous_solutions());

    // Prepare previous solutions transfer
    previous_solutions_transfer.reserve(previous_solutions.size());
    for (unsigned int i = 0; i < previous_solutions.size(); ++i)
      {
        previous_solutions_transfer.emplace_back(
          parallel::distributed::
            SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>(
              this->dof_handler));
      }

    // Check the value of interface sharpness
    if (simulation_parameters.multiphysics.vof_parameters.sharpening
          .interface_sharpness < 1.0)
      this->pcout
        << "Warning: interface sharpness values smaller than 1 smooth the interface instead of sharpening it."
        << std::endl
        << "The interface sharpness value should be set between 1 and 2"
        << std::endl;


    // Change the behavior of the timer for situations when you don't want
    // outputs
    if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
      this->computing_timer.disable_output();
  }

  /**
   * @brief VOF - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  ~VolumeOfFluid()
  {}


  /**
   * @brief Call for the assembly of the matrix and the right-hand side.
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver
   */
  void
  assemble_matrix_and_rhs();

  /**
   * @brief Call for the assembly of the right-hand side
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver
   */
  void
  assemble_rhs();

  /**
   * @brief Attach the solution vector to the DataOut provided. This function
   * enable the auxiliary physics to output their solution via the core solver.
   */
  void
  attach_solution_to_output(DataOut<dim> &data_out) override;

  /**
   * @brief Calculates the L2 error of the solution
   */
  double
  calculate_L2_error();

  /**
   * @brief Calculates the volume and mass for a given fluid phase.
   * Used for conservation monitoring.
   *
   * @param solution VOF solution (phase fraction)
   *
   * @param current_solution_fd current solution for the fluid dynamics
   *
   * @param monitored_fluid Fluid indicator (fluid0 or fluid1) corresponding to
   * the phase of interest.
   */
  template <typename VectorType>
  void
  calculate_volume_and_mass(const TrilinosWrappers::MPI::Vector &solution,
                            const VectorType &current_solution_fd,
                            const Parameters::FluidIndicator monitored_fluid);

  /**
   * @brief Calculates the barycenter of the fluid and its velocity
   *
   * @param solution VOF solution
   *
   * @param solution Fluid dynamics solution
   *
   */
  template <typename VectorType>
  std::pair<Tensor<1, dim>, Tensor<1, dim>>
  calculate_barycenter(const TrilinosWrappers::MPI::Vector &solution,
                       const VectorType &current_solution_fd);


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
   * @brief Compute the Kelly error estimator on the phase parameter for mesh refinement.
   * See :
   * https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html
   * for more information on the Kelly error estimator.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  void
  compute_kelly(const std::pair<const Parameters::MeshAdaptation::Variable,
                                Parameters::MultipleAdaptationParameters> &ivar,
                dealii::Vector<float> &estimated_error_per_cell) override;

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
   * only support imposing nodal values, but some physics additionnaly support
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
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated or not
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   * NB : dof_handler and present_solution are passed to the multiphysics
   * interface at the end of the setup_dofs method
   */
  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return dof_handler;
  }
  TrilinosWrappers::MPI::Vector &
  get_evaluation_point() override
  {
    return evaluation_point;
  }
  TrilinosWrappers::MPI::Vector &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  }
  TrilinosWrappers::MPI::Vector &
  get_newton_update() override
  {
    return newton_update;
  }
  TrilinosWrappers::MPI::Vector &
  get_present_solution() override
  {
    return present_solution;
  }
  TrilinosWrappers::MPI::Vector &
  get_system_rhs() override
  {
    return system_rhs;
  }
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }
  DoFHandler<dim> *
  get_projected_phase_fraction_gradient_dof_handler()
  {
    return &projected_phase_fraction_gradient_dof_handler;
  }
  DoFHandler<dim> *
  get_curvature_dof_handler()
  {
    return &curvature_dof_handler;
  }
  TrilinosWrappers::MPI::Vector *
  get_projected_phase_fraction_gradient_solution()
  {
    return &present_projected_phase_fraction_gradient_solution;
  }
  TrilinosWrappers::MPI::Vector *
  get_curvature_solution()
  {
    return &present_curvature_solution;
  }

private:
  /**
   *  @brief Assembles the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;


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
   * See the documentation for VOFScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    VOFScratchData<dim>                                  &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for VOFScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    VOFScratchData<dim>                                  &scratch_data,
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
   * @brief Limit the phase fractions between 0 and 1. This is necessary before interface sharpening.
   * More information can be found in step_41 of deal.II tutorials:
   * https://www.dealii.org/current/doxygen/deal.II/step_41.html
   */
  void
  update_solution_and_constraints(TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Assemble the system for interface sharpening
   *  * This function assembles the weak form of:
   *  $$ \Phi = c ^ (1 - \alpha) * (\phi ^ \alpha)  if 0 <=
   * \phi <= c $$
   *  $$ \Phi = 1 - (1 - c) ^ (1 - \alpha) * (1 - \phi) ^ \alpha  if c <
   * \phi <= 1 $$
   * Reference for sharpening method
   * https://www.sciencedirect.com/science/article/pii/S0045782500002000
   *
   * @param solution VOF solution (phase fraction)
   *
   * @param sharpening_threshold Interface sharpening threshold that represents
   * the mass conservation level
   *
   */
  void
  assemble_L2_projection_interface_sharpening(
    TrilinosWrappers::MPI::Vector &solution,
    const double                   sharpening_threshold);

  /**
   * @brief Solves the assembled system to sharpen the interface. The linear_solver_tolerance
   * is hardcoded = 1e-15, and an ILU Preconditioner is used. After solving the
   * system, this function overwrites the solution with the sharpened solution
   *
   * @param solution VOF solution (phase fraction)
   */
  void
  solve_interface_sharpening(TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Assembles a mass_matrix which is used in update_solution_and_constraints function
   * to limit the phase fractions of cells in the range of [0,1] before
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
   * monitored, as given in the prm (VOF, subsection monitoring), by binary
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
   * @param solution VOF solution (phase fraction)
   *
   * @param sharpening_threshold Interface sharpening threshold that represents the mass conservation level
   *
   * @param sharpen_previous_solutions Boolean true if sharpening is applied on the present
   * and past solutions, false if the sharpening is applied on the given
   * solution vector only. Used to determine the sharpening threshold by binary
   * search.
   */
  void
  sharpen_interface(TrilinosWrappers::MPI::Vector &solution,
                    const double                   sharpening_threshold,
                    const bool                     sharpen_previous_solutions);

  /**
   * @brief Carries out the smoothing phase fraction with a projection step (to avoid a staircase interface).
   */
  void
  smooth_phase_fraction();

  /**
   * @brief Carries out finding the gradients of phase fraction. Obtained gradients of phase
   * fraction is used in find_projected_interface_curvature to find interface
   * curvature (k).
   */
  void
  find_projected_phase_fraction_gradient();

  /**
   * @brief Carries out finding the interface curvature.
   */
  void
  find_projected_interface_curvature();

  /**
   * @brief Assembles the matrix and rhs for calculation of a smooth phase fraction using a projection.
   *
   * @param solution VOF solution (phase fraction)
   */
  void
  assemble_projection_phase_fraction(TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Solves smooth phase fraction system.
   * @param solution VOF solution (phase fraction)
   */
  void
  solve_projection_phase_fraction(TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Assembles the matrix and rhs for calculation of projected phase fraction gradient (pfg).
   *
   * Solves:
   * $$ v . \psi + \eta * \nabla v . \nabla \psi = v . \nabla \phi $$
   * where $$v$$, $$\psi$$, $$\eta$$, and $$\phi$$ are test function, fpg,
   *
   * value, and phase fraction.
   *
   * @param solution VOF solution (phase fraction)
   */
  void
  assemble_projected_phase_fraction_gradient_matrix_and_rhs(
    TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Solves phase fraction gradient system.
   */
  void
  solve_projected_phase_fraction_gradient();

  /**
   * @brief Assembles the matrix and rhs for calculation of the curvature.
   *
   * Solves:
   * $$ v * k + \eta * \nabla v . \nabla k = \nabla v . (\psi / |\psi|) $$
   * where $$v$$, $$psi$$, $$eta$$, and $$k$$ are test function, fpg, filter
   * value, and curvature.
   *
   * @param present_projected_phase_fraction_gradient_solution
   */
  void
  assemble_curvature_matrix_and_rhs(
    TrilinosWrappers::MPI::Vector
      &present_projected_phase_fraction_gradient_solution);

  /**
   * @brief Solves curvature system.
   */
  void
  solve_curvature();

  /**
   * @brief Applies filter on phase fraction values.
   */
  void
  apply_phase_filter();


  TrilinosWrappers::MPI::Vector nodal_phase_fraction_owned;

  MultiphysicsInterface<dim> *multiphysics;

  TimerOutput computing_timer;

  const SimulationParameters<dim> &simulation_parameters;

  // Core elements for the VOF simulation
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl> simulation_control;
  DoFHandler<dim>                    dof_handler;
  DoFHandler<dim> projected_phase_fraction_gradient_dof_handler;
  DoFHandler<dim> curvature_dof_handler;
  std::shared_ptr<FiniteElement<dim>> fe;
  std::shared_ptr<FESystem<dim>>      fe_projected_phase_fraction_gradient;
  std::shared_ptr<FiniteElement<dim>> fe_curvature;
  ConvergenceTable                    error_table;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>    mapping;
  std::shared_ptr<Mapping<dim>>    projected_phase_fraction_gradient_mapping;
  std::shared_ptr<Mapping<dim>>    curvature_mapping;
  std::shared_ptr<Quadrature<dim>> cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  // Solution storage
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  TrilinosWrappers::MPI::Vector  evaluation_point;
  TrilinosWrappers::MPI::Vector  local_evaluation_point;
  TrilinosWrappers::MPI::Vector  newton_update;
  TrilinosWrappers::MPI::Vector  present_solution;
  TrilinosWrappers::MPI::Vector  system_rhs;
  AffineConstraints<double>      nonzero_constraints;
  AffineConstraints<double>      bounding_constraints;
  AffineConstraints<double>      zero_constraints;
  TrilinosWrappers::SparseMatrix system_matrix;
  TrilinosWrappers::MPI::Vector  solution_pw;
  TrilinosWrappers::MPI::Vector  filtered_solution;

  // Previous solutions vectors
  std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;

  // Solution transfer classes
  std::shared_ptr<
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>
    solution_transfer;
  std::vector<
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>
    previous_solutions_transfer;

  // Phase fraction matrices for interface sharpening
  TrilinosWrappers::SparseMatrix system_matrix_phase_fraction;
  TrilinosWrappers::SparseMatrix complete_system_matrix_phase_fraction;
  TrilinosWrappers::MPI::Vector  system_rhs_phase_fraction;
  TrilinosWrappers::MPI::Vector  complete_system_rhs_phase_fraction;
  TrilinosWrappers::SparseMatrix mass_matrix_phase_fraction;

  // Projected phase fraction gradient (pfg) solution
  TrilinosWrappers::MPI::Vector
           present_projected_phase_fraction_gradient_solution;
  IndexSet locally_owned_dofs_projected_phase_fraction_gradient;
  IndexSet locally_relevant_dofs_projected_phase_fraction_gradient;
  AffineConstraints<double> projected_phase_fraction_gradient_constraints;
  TrilinosWrappers::MPI::Vector
    nodal_projected_phase_fraction_gradient_relevant;
  TrilinosWrappers::MPI::Vector nodal_projected_phase_fraction_gradient_owned;

  TrilinosWrappers::SparseMatrix
                                system_matrix_projected_phase_fraction_gradient;
  TrilinosWrappers::MPI::Vector system_rhs_projected_phase_fraction_gradient;

  // Projected curvature solution
  TrilinosWrappers::MPI::Vector present_curvature_solution;
  IndexSet                      locally_owned_dofs_curvature;
  IndexSet                      locally_relevant_dofs_curvature;
  AffineConstraints<double>     curvature_constraints;
  TrilinosWrappers::MPI::Vector nodal_curvature_relevant;
  TrilinosWrappers::MPI::Vector nodal_curvature_owned;

  std::vector<Tensor<1, dim>> projected_phase_fraction_gradient_values;
  std::vector<double>         curvature_values;

  TrilinosWrappers::SparseMatrix                     system_matrix_curvature;
  TrilinosWrappers::MPI::Vector                      system_rhs_curvature;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Lower and upper bounds of phase fraction
  const double phase_upper_bound = 1.0;
  const double phase_lower_bound = 0.0;

  // Conservation Analysis
  TableHandler table_monitoring_vof;
  double       volume_monitored;
  double       mass_monitored;
  double       mass_first_iteration;
  double       sharpening_threshold;

  // Barycenter analysis
  TableHandler table_barycenter;

  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<VOFAssemblerBase<dim>>> assemblers;

  // Phase fraction filter
  std::shared_ptr<VolumeOfFluidFilterBase> filter;
};



#endif
