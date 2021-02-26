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
 * Implementation of heat transfer as an auxiliary physics.
 * This heat equation is weakly coupled to the velocity field.
 * Equation solved:
 * rho * Cp * (dT/dt + u.gradT) = k div(gradT) + nu/rho * (gradu : gradu)
 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020-
 */

#ifndef lethe_heat_transfer_h
#define lethe_heat_transfer_h

#include <deal.II/base/convergence_table.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <core/simulation_control.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/base/quadrature_lib.h>
#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>


template <int dim>
class HeatTransfer : public AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>
{
public:
  HeatTransfer<dim>(MultiphysicsInterface<dim> *     multiphysics_interface,
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
    , solution_transfer_m1(dof_handler)
    , solution_transfer_m2(dof_handler)
    , solution_transfer_m3(dof_handler)
  {
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe = std::make_shared<Simplex::FE_P<dim>>(
          simulation_parameters.fem_parameters.temperature_order);
        temperature_mapping = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature =
          std::make_shared<Simplex::QGauss<dim>>(fe->degree + 1);
        face_quadrature =
          std::make_shared<Simplex::QGauss<dim - 1>>(fe->degree + 1);
        error_quadrature =
          std::make_shared<Simplex::QGauss<dim>>(fe->degree + 2);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe = std::make_shared<FE_Q<dim>>(
          simulation_parameters.fem_parameters.temperature_order);
        temperature_mapping = std::make_shared<MappingQ<dim>>(
          fe->degree, simulation_parameters.fem_parameters.qmapping_all);
        cell_quadrature  = std::make_shared<QGauss<dim>>(fe->degree + 1);
        face_quadrature  = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
        error_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 2);
      }
  }

  /**
   * @brief Call for the assembly of the matrix and the right-hand side.
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  /**
   * @brief Call for the assembly of the right-hand side
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  /**
   * @brief Attach the solution vector to the DataOut provided. This function
   * enable the auxiliary physics to output their solution via the core solver.
   */
  virtual void
  attach_solution_to_output(DataOut<dim> &data_out) override;


  /**
   * @brief Calculates the L2 error of the solution
   */
  double
  calculate_L2_error();


  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation() override;

  /**
   * @brief Carry out the operations require to finish a time step correctly. This
   * includes setting the previous values
   */
  virtual void
  finish_time_step() override;

  /**
   * @brief Rearrange vector solution correctly for transient simulations
   */
  virtual void
  percolate_time_vectors() override;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess(bool first_iteration) override;


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  virtual void
  pre_mesh_adaptation();

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  virtual void
  post_mesh_adaptation();

  /**
   * @brief Prepares Heat Transfer to write checkpoint
   */
  virtual void
  write_checkpoint() override;

  /**
   * @brief Allows Heat Transfer to set-up solution vector from checkpoint file;
   */
  virtual void
  read_checkpoint() override;


  /**
   * @brief Returns the dof_handler of the heat transfer physics
   */
  virtual const DoFHandler<dim> &
  get_dof_handler() override
  {
    return dof_handler;
  };

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
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
  virtual void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true);


  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   *
   * @param number_physic_current Indicates the number associated with the physic currently solved
   * default value = 0, meaning only one physic, generally associated with the
   * flow equations, is solved by default
   */
  virtual TrilinosWrappers::MPI::Vector &
  get_evaluation_point() override
  {
    return evaluation_point;
  }
  virtual TrilinosWrappers::MPI::Vector &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  }
  virtual TrilinosWrappers::MPI::Vector &
  get_newton_update() override
  {
    return newton_update;
  }
  virtual TrilinosWrappers::MPI::Vector &
  get_present_solution() override
  {
    return present_solution;
  }
  virtual TrilinosWrappers::MPI::Vector &
  get_system_rhs() override
  {
    return system_rhs;
  }
  virtual AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }


private:
  template <bool assemble_matrix>
  void
  assemble_system(const Parameters::SimulationControl::TimeSteppingMethod
                    time_stepping_method);

  MultiphysicsInterface<dim> *     multiphysics;
  const SimulationParameters<dim> &simulation_parameters;


  // Core elements for the heat transfer simulation
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl> simulation_control;
  DoFHandler<dim>                    dof_handler;

  std::shared_ptr<FiniteElement<dim>> fe;
  ConvergenceTable                    error_table;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>        temperature_mapping;
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


  // Past solution vectors
  TrilinosWrappers::MPI::Vector solution_m1;
  TrilinosWrappers::MPI::Vector solution_m2;
  TrilinosWrappers::MPI::Vector solution_m3;

  // Solution transfer classes
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    solution_transfer;
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    solution_transfer_m1;
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    solution_transfer_m2;
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    solution_transfer_m3;

  // Reference for GGLS https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.2324
  // Warning, this GGLS implementation is valid only for Linear elements
  // Quad elements will be lacking the third derivative of the diffusion
  // operator Whether this affects or not the final result is unclear to me at
  // the moment. Additionnaly, this formulation does not use the gradient of the
  // source term. The same applies, I have no clue if this is detrimental or not
  // to the solution since anyway the GGLS term scales as h^(order+1)
  const bool GGLS = true;
};


#endif
