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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef lethe_navier_stokes_base_h
#define lethe_navier_stokes_base_h

// Dealii Includes

// Base
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// Lac
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

// Grid
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Numerics
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Distributed
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>


// Lethe Includes
#include <core/bdf.h>
#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/newton_non_linear_solver.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>
#include <solvers/flow_control.h>
#include <solvers/postprocessing_cfd.h>
#include <solvers/postprocessing_velocities.h>

#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"

// Std
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * A base class for all the Navier-Stokes equation
 * This class regroups common facilities that are shared by all
 * the Navier-Stokes implementations to reduce code multiplicity
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @tparam VectorType  The Vector type used for the solvers
 *
 * @tparam DofsType the type of dof storage indices
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim, typename VectorType, typename DofsType>
class NavierStokesBase : public PhysicsSolver<VectorType>
{
protected:
  NavierStokesBase(NavierStokesSolverParameters<dim> &nsparam);

  virtual ~NavierStokesBase()
  {}

  /**
   * @brief calculate_forces
   * Post-processing function
   * Calculate forces acting on each boundary condition
   */
  void
  postprocessing_forces(const VectorType &evaluation_point);

  /**
   * @brief calculate_L2_error
   * @return L2 norm of the error for velocity and pressure
   * Post-processing function
   * Calculates the L2 norm of the error
   */
  std::pair<double, double>
  calculate_L2_error(const VectorType &evaluation_point);

  /**
   * @brief calculate_torques
   * Post-processing function
   * Calculate torque acting on each boundary condition
   */
  void
  postprocessing_torques(const VectorType &evaluation_point);

  /**
   * @brief calculate_flow_rate
   * @return inlet flow rate and area
   * Post-processing function
   * Calculate the volumetric flow rate on the inlet boundary
   */
  void
  postprocessing_flow_rate(const VectorType &evaluation_point);

  /**
   * @brief dynamic_flow_control
   * If set to enable, dynamic_flow_control allows to control the flow by
   * calculate a beta coefficient at each time step added to the force of the
   * source term for gls_navier_stokes solver.
   */
  void
  dynamic_flow_control();

  /**
   * @brief finish_time_step
   * Finishes the time step
   * Post-processing and time stepping
   */
  virtual void
  finish_time_step();

  /**
   * @brief finish_simulation
   * Finishes the simulation by calling all
   * the post-processing elements that are required
   */
  void
  finish_simulation();

  /**
   * @brief iterate
   * Do a regular CFD iteration
   */
  virtual void
  iterate();


  /**
   * @brief First iteration
   * Do the first CFD iteration
   */
  virtual void
  first_iteration();

  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) = 0;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;

  void
  refine_mesh();

  void
  refine_mesh_kelly();

  void
  refine_mesh_uniform();

  /**
   * @brief postprocess
   * Post-process after an iteration
   */
  virtual void
  postprocess(bool firstIter);

  /**
   * @brief read_checkpoint
   */
  virtual void
  read_checkpoint();

  /**
   * @brief set_nodal_values
   * Set the nodal values of velocity and pressure
   */
  void
  set_nodal_values();

  /**
   * @brief set_periodicity
   *
   * Initialize the periodic boundary conditions
   */
  virtual void
  setup_dofs() = 0;

  /**
   * @brief write_checkpoint
   */
  virtual void
  write_checkpoint();

  /**
   * @brief write_output_results
   * Post-processing as parallel VTU files
   */
  void
  write_output_results(const VectorType &solution);

  /**
   * @brief write_output_forces
   * Writes the forces per boundary condition to a text file output
   */
  virtual void
  output_field_hook(DataOut<dim> &);

  void
  write_output_forces();

  /**
   * @brief write_output_torques
   * Writes the torques per boundary condition to a text file output
   */
  void
  write_output_torques();

  // Member variables
  DofsType locally_owned_dofs;
  DofsType locally_relevant_dofs;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  DoFHandler<dim>                                              dof_handler;
  FESystem<dim>                                                fe;

  TimerOutput computing_timer;

  NavierStokesSolverParameters<dim> nsparam;
  PVDHandler                        pvdhandler;

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

  // Inlet flow rate and area
  std::pair<double, double> flow_rate;

  // Dynamic flow control
  FlowControl<dim> flow_control;
  Tensor<1, dim>   beta;

  // Constraints for Dirichlet boundary conditions
  AffineConstraints<double> zero_constraints;

  // Solution vectors
  VectorType solution_m1;
  VectorType solution_m2;
  VectorType solution_m3;

  // Finite element order used
  const unsigned int velocity_fem_degree;
  const unsigned int pressure_fem_degree;
  unsigned int       number_quadrature_points;

  // Simulation control for time stepping and I/Os
  std::shared_ptr<SimulationControl> simulation_control;
  // SimulationControl simulationControl;

  // Post-processing variables
  TableHandler                                 enstrophy_table;
  TableHandler                                 kinetic_energy_table;
  AverageVelocities<dim, VectorType, DofsType> average_velocities;
  VectorType                                   average_solution;

  // Convergence Analysis
  ConvergenceTable error_table;

  // Force analysis
  std::vector<Tensor<1, dim>> forces_on_boundaries;
  std::vector<Tensor<1, 3>>   torques_on_boundaries;
  std::vector<TableHandler>   forces_tables;
  std::vector<TableHandler>   torques_tables;
};

#endif
