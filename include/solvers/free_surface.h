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
 * Implementation of free surface with a Volume of Fluid method.
 * Two fluid formulation. The phase indicator "phase" is equal to 0
 * in one fluid and 1 in the other. The free surface is located
 * where "phase" is equal to 0.5.
 *
 * Author: Jeanne Joachim, Polytechnique Montreal, 2021
 */

#ifndef lethe_free_surface_h
#define lethe_free_surface_h

#include <core/bdf.h>
#include <core/simulation_control.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>

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

#include <deal.II/numerics/error_estimator.h>


template <int dim>
class FreeSurface : public AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>
{
public:
  /**
   * @brief FreeSurface - Base constructor.
   */
  FreeSurface<dim>(MultiphysicsInterface<dim> *     multiphysics_interface,
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
  }

  /**
   * @brief FreeSurface - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  ~FreeSurface()
  {}

  /**
   * @brief Call for the assembly of the matrix
   */
  virtual void
  assemble_system_matrix()
  {
    assemble_matrix_and_rhs();
  }

  /**
   * @brief Call for the assembly of the right-hand side
   */
  virtual void
  assemble_system_rhs()
  {
    assemble_rhs();
  }


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
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  void
  finish_simulation();

  /**
   * @brief Carry out the operations required to finish a time step correctly.
   */
  void
  finish_time_step();

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  void
  percolate_time_vectors();

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  void
  postprocess(bool first_iteration);


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
  write_checkpoint();


  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  void
  read_checkpoint();

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  void
  setup_dofs();

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  void
  set_initial_conditions();

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

private:
  template <bool assemble_matrix>
  void
  assemble_system();

  MultiphysicsInterface<dim> *     multiphysics;
  const SimulationParameters<dim> &simulation_parameters;


  // Core elements for the free surface simulation
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl> simulation_control;
  DoFHandler<dim>                    dof_handler;

  std::shared_ptr<FiniteElement<dim>> fe;
  ConvergenceTable                    error_table;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>        mapping;
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
  AffineConstraints<double>      zero_constraints;
  TrilinosWrappers::SparseMatrix system_matrix;


  // Previous solutions vectors
  std::vector<TrilinosWrappers::MPI::Vector> previous_solutions;

  // Solution transfer classes
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    solution_transfer;
  std::vector<
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>
    previous_solutions_transfer;

  // Enable DCDD shock capturing scheme
  const bool DCDD = true;
};



#endif
