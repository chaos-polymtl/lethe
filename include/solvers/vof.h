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

template <int dim>
class VolumeOfFluid
  : public AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>
{
public:
  /**
   * @brief VOF - Base constructor.
   */
  VolumeOfFluid<dim>(
    MultiphysicsInterface<dim> *     multiphysics_interface,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control)
    : AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>(
        p_simulation_parameters.non_linear_solver)
    , multiphysics(multiphysics_interface)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
    , pfg_dof_handler(*triangulation)
    , curvature_dof_handler(*triangulation)
    , solution_transfer(dof_handler)
  {
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe              = std::make_shared<FE_SimplexP<dim>>(1);
        mapping         = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
        face_quadrature =
          std::make_shared<QGaussSimplex<dim - 1>>(fe->degree + 1);
        error_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 2);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe      = std::make_shared<FE_Q<dim>>(1);
        mapping = std::make_shared<MappingQ<dim>>(
          fe->degree, simulation_parameters.fem_parameters.qmapping_all);
        fe_pfg = std::make_shared<FESystem<dim>>(FE_Q<dim>(fe->degree), dim);
        fe_curvature = std::make_shared<FE_Q<dim>>(1);
        pfg_mapping  = std::make_shared<MappingQ<dim>>(
          fe_pfg->degree, simulation_parameters.fem_parameters.qmapping_all);
        curvature_mapping = std::make_shared<MappingQ<dim>>(
          fe_curvature->degree,
          simulation_parameters.fem_parameters.qmapping_all);
        cell_quadrature  = std::make_shared<QGauss<dim>>(fe->degree + 1);
        face_quadrature  = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
        error_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 2);
      }

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
    if (simulation_parameters.interface_sharpening.interface_sharpness < 1.0)
      this->pcout
        << "Warning: interface sharpness values smaller than 1 smooth the interface instead of sharpening it."
        << std::endl
        << "The interface sharpness value should be set between 1 and 2"
        << std::endl;
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
  attach_solution_to_output(DataOut<dim> &data_out);

  /**
   * @brief Calculates the L2 error of the solution
   */
  double
  calculate_L2_error();

  /**
   * @brief Calculates the volume for the fluid phase with given id_fluid_monitored.
   * Used for conservation monitoring.
   */
  double
  calculate_volume(int id_fluid_monitored);

  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  void
  finish_simulation() override;

  /**
   * @brief Carry out the operations required to finish a time step correctly.
   */
  void
  finish_time_step() override;

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Carry out modifications on the auxiliary physic solution.
   * Used in vof method for interface sharpening and wetting/peeling.
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
  pre_mesh_adaptation();

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  void
  post_mesh_adaptation();

  /**
   * @brief Compute the Kelly error estimator on the phase parameter for mesh refinement.
   * See :
   * https://www.dealii.org/current/doxygen/deal.II/classKellyErrorEstimator.html
   * for more information on the Kelly error estimator.
   */
  void
  compute_kelly(dealii::Vector<float> &estimated_error_per_cell);

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
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  void
  set_initial_conditions() override;

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
                      const bool renewed_matrix = true);

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

  // enum class used for peeling/wetting
  enum class PhaseChange
  {
    peeling,
    wetting
  };

private:
  /**
   *  @brief Assembles the matrix associated with the solver
   */
  void
  assemble_system_matrix();

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs();


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
    VOFScratchData<dim> &                                 scratch_data,
    StabilizedMethodsCopyData &                           copy_data);

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
    VOFScratchData<dim> &                                 scratch_data,
    StabilizedMethodsCopyData &                           copy_data);

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
   */
  void
  assemble_L2_projection_interface_sharpening(
    TrilinosWrappers::MPI::Vector &solution);

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
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  /**
   * @brief Modification of the solution
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param i_bc peeling-wetting boundary index
   *
   * @param current_solution_fd current solution for the fluid dynamics
   */
  template <typename VectorType>
  void
  apply_peeling_wetting(const unsigned int i_bc,
                        const VectorType & current_solution_fd);

  /**
   * @brief Change cell phase, small method called to avoid code repetition and reduce sloppy
   * error likelihood in apply_peeling_wetting.
   *
   * @param type a parameter of class PhaseChange (see below) stating the needed change
   *
   * @param new_phase the new phase value for the cell (0 or 1)
   *
   * @param solution_pw VOF solution after peeling and wetting corrections are applied
   *
   * @param dof_indices_vof local index for the VOF solution
   */
  void
  change_cell_phase(
    const PhaseChange &                         type,
    const unsigned int &                        new_phase,
    TrilinosWrappers::MPI::Vector &             solution_pw,
    const std::vector<types::global_dof_index> &dof_indices_vof);

  /**
   * @brief Carries out interface sharpening. It is called in the modify solution function.
   */
  void
  sharpen_interface();

  /**
   * @brief Carries out finding the gradients of phase fraction. Obtained gradients of phase
   * fraction is used in find_filtered_interface_curvature to find interface
   * curvature (k).
   */
  void
  find_filtered_pfg();

  /**
   * @brief Carries out finding the interface curvature.
   */
  void
  find_filtered_interface_curvature();

  /**
   * @brief Assembles the matrix and rhs for calculation of filtered phase gradient (fpg).
   *
   * Solves:
   * $$ v . \psi + \eta * \nabla v . \nabla \psi = v . \nabla \phi $$
   * where $$v$$, $$\psi$$, $$\eta$$, and $$\phi$$ are test function, fpg,
   * filter value, and phase fraction.
   *
   * @param solution VOF solution (phase fraction)
   */
  void
  assemble_pfg_matrix_and_rhs(TrilinosWrappers::MPI::Vector &solution);

  /**
   * @brief Solves phase fraction gradient system.
   */
  void
  solve_pfg();

  /**
   * @brief Assembles the matrix and rhs for calculation of the curvature.
   *
   * Solves:
   * $$ v * k + \eta * \nabla v . \nabla k = \nabla v . (\psi / |\psi|) $$
   * where $$v$$, $$psi$$, $$eta$$, and $$k$$ are test function, fpg, filter
   * value, and curvature.
   *
   * @param present_pfg_solution
   */
  void
  assemble_curvature_matrix_and_rhs(
    TrilinosWrappers::MPI::Vector &present_pfg_solution);

  /**
   * @brief Solves curvature system.
   */
  void
  solve_curvature();

  TrilinosWrappers::MPI::Vector nodal_phase_fraction_owned;

  MultiphysicsInterface<dim> *     multiphysics;
  const SimulationParameters<dim> &simulation_parameters;

  // Core elements for the VOF simulation
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl>  simulation_control;
  DoFHandler<dim>                     dof_handler;
  DoFHandler<dim>                     pfg_dof_handler;
  DoFHandler<dim>                     curvature_dof_handler;
  std::shared_ptr<FiniteElement<dim>> fe;
  std::shared_ptr<FESystem<dim>>      fe_pfg;
  std::shared_ptr<FiniteElement<dim>> fe_curvature;
  ConvergenceTable                    error_table;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Mapping<dim>>        pfg_mapping;
  std::shared_ptr<Mapping<dim>>        curvature_mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;
  std::shared_ptr<Quadrature<dim>>     error_quadrature;

  // Solution storage:
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

  // Previous solutions vectors
  std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;
  std::vector<TrilinosWrappers::MPI::Vector> solution_stages;

  // Solution transfer classes
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
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

  // Peeling/Wetting analysis
  TrilinosWrappers::MPI::Vector marker_pw;

  // Filtered phase fraction gradient (pfg) solution
  TrilinosWrappers::MPI::Vector present_pfg_solution;
  IndexSet                      locally_owned_dofs_pfg;
  IndexSet                      locally_relevant_dofs_pfg;
  AffineConstraints<double>     pfg_constraints;
  TrilinosWrappers::MPI::Vector nodal_pfg_relevant;
  TrilinosWrappers::MPI::Vector nodal_pfg_owned;

  TrilinosWrappers::SparseMatrix system_matrix_pfg;
  TrilinosWrappers::MPI::Vector  system_rhs_pfg;

  // Filtered curvature solution
  TrilinosWrappers::MPI::Vector present_curvature_solution;
  IndexSet                      locally_owned_dofs_curvature;
  IndexSet                      locally_relevant_dofs_curvature;
  AffineConstraints<double>     curvature_constraints;
  TrilinosWrappers::MPI::Vector nodal_curvature_relevant;
  TrilinosWrappers::MPI::Vector nodal_curvature_owned;

  std::vector<Tensor<1, dim>> pfg_values;
  std::vector<double>         curvature_values;

  TrilinosWrappers::SparseMatrix system_matrix_curvature;
  TrilinosWrappers::MPI::Vector  system_rhs_curvature;

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Lower and upper bounds of phase fraction
  const double phase_upper_bound = 1.0;
  const double phase_lower_bound = 0.0;

  // Conservation Analysis
  TableHandler table_monitoring_vof;

  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<VOFAssemblerBase<dim>>> assemblers;
};



#endif
