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

 *
 * Implementation of tracer as an auxiliary physics.
 * Equation solved:
 * dT/dt +  u * gradT = D * div(grad T) + f
 * with T the tracer function, D the diffusivity and f the forcing
 *
 */

#ifndef lethe_tracer_h
#define lethe_tracer_h

#include <core/bdf.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/tracer_assemblers.h>
#include <solvers/tracer_scratch_data.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>


template <int dim>
class Tracer : public AuxiliaryPhysics<dim, GlobalVectorType>
{
public:
  Tracer(MultiphysicsInterface<dim>      *multiphysics_interface,
         const SimulationParameters<dim> &p_simulation_parameters,
         std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                            p_triangulation,
         std::shared_ptr<SimulationControl> p_simulation_control)
    : AuxiliaryPhysics<dim, GlobalVectorType>(
        p_simulation_parameters.non_linear_solver.at(PhysicsID::tracer))
    , multiphysics(multiphysics_interface)
    , computing_timer(p_triangulation->get_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
  {
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe = std::make_shared<FE_SimplexP<dim>>(
          simulation_parameters.fem_parameters.tracer_order);
        mapping         = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe = std::make_shared<FE_Q<dim>>(
          simulation_parameters.fem_parameters.tracer_order);
        mapping         = std::make_shared<MappingQ<dim>>(fe->degree);
        cell_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
      }

    // Allocate solution transfer
    solution_transfer = std::make_shared<
      parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>(
      dof_handler);

    // Set size of previous solutions using BDF schemes information
    previous_solutions.resize(maximum_number_of_previous_solutions());

    // Prepare previous solutions transfer
    previous_solutions_transfer.reserve(previous_solutions.size());
    for (unsigned int i = 0; i < previous_solutions.size(); ++i)
      {
        previous_solutions_transfer.emplace_back(
          parallel::distributed::SolutionTransfer<dim, GlobalVectorType>(
            this->dof_handler));
      }

    // Change the behavior of the timer for situations when you don't want
    // outputs
    if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
      this->computing_timer.disable_output();
  }

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
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  void
  finish_simulation() override;

  /**
   * @brief Rearrange vector solution correctly for transient simulations
   */
  void
  percolate_time_vectors() override;

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
   * @brief Compute the Kelly error estimator for mesh refinement.
   * NB : not implemented for the tracer parameter for now.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  void
  compute_kelly(
    const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
      & /*ivar*/,
    dealii::Vector<float> & /*estimated_error_per_cell*/) override
  {
    return;
  }

  /**
   * @brief Prepares Heat Transfer to write checkpoint
   */
  void
  write_checkpoint() override;

  /**
   * @brief Allows tracer physics to set-up solution vector from checkpoint file;
   */
  void
  read_checkpoint() override;


  /**
   * @brief Returns the dof_handler of the tracer physics
   */
  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return dof_handler;
  }

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
   * @brief Update non zero constraints if the boundary is time-dependent
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
    return present_solution;
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


  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  void
  output_newton_update_norms(const unsigned int display_precision) override
  {
    this->pcout << std::setprecision(display_precision)
                << "\t||dx||_L2 = " << std::setw(6) << newton_update.l2_norm()
                << std::setw(6)
                << "\t||dx||_Linfty = " << std::setprecision(display_precision)
                << newton_update.linfty_norm() << std::endl;
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
   * See the documentation for TracerScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    TracerScratchData<dim>                               &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for TracerScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    TracerScratchData<dim>                               &scratch_data,
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
   * @brief Calculate tracer statistics : Max, min, average and standard-deviation
   */
  void
  calculate_tracer_statistics();

  /**
   * @brief Writes the tracer statistics to an output file
   */
  void
  write_tracer_statistics();

  MultiphysicsInterface<dim> *multiphysics;

  TimerOutput computing_timer;

  const SimulationParameters<dim> &simulation_parameters;


  // Core elements for the tracer
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl> simulation_control;
  DoFHandler<dim>                    dof_handler;

  // Finite element spce
  std::shared_ptr<FiniteElement<dim>> fe;
  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>    mapping;
  std::shared_ptr<Quadrature<dim>> cell_quadrature;


  ConvergenceTable error_table;

  // Solution storage:
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  GlobalVectorType               evaluation_point;
  GlobalVectorType               local_evaluation_point;
  GlobalVectorType               newton_update;
  GlobalVectorType               present_solution;
  GlobalVectorType               system_rhs;
  AffineConstraints<double>      nonzero_constraints;
  AffineConstraints<double>      zero_constraints;
  TrilinosWrappers::SparseMatrix system_matrix;


  // Previous solutions vectors
  std::vector<GlobalVectorType> previous_solutions;

  // Solution transfer classes
  std::shared_ptr<
    parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>
    solution_transfer;
  std::vector<parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>
    previous_solutions_transfer;

  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<TracerAssemblerBase<dim>>> assemblers;

  // Tracer statistics table
  TableHandler statistics_table;
};


#endif
