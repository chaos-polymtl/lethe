﻿/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#ifndef LETHE_GDNS_H
#define LETHE_GDNS_H

#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#include "navier_stokes_base.h"

using namespace dealii;

/**
 * A solver class for the Steady-Sate  Navier-Stokes equation using Grad-Div
 * stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim>
class GDNavierStokesSolver
  : public NavierStokesBase<dim,
                            TrilinosWrappers::MPI::BlockVector,
                            std::vector<IndexSet>>
{
public:
  GDNavierStokesSolver(NavierStokesSolverParameters<dim> &nsparam,
                       const unsigned int                 degreeVelocity,
                       const unsigned int                 degreePressure);
  ~GDNavierStokesSolver();

  void
  solve();

private:
  void
  assemble_matrix_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) override;

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme>
  void
  assembleGD();
  void
  assemble_L2_projection();

  void
  set_nodal_values();

  virtual void
  setup_dofs();

  void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool                             restart = false);
  void
  set_solution_vector(double value);

  /**
   * Solver for the L2 Projection linear system
   */

  void
  solve_L2_system(bool   initial_step,
                  double relative_residual,
                  double minimum_residual);

  /**
   * Solver for the NS linear system of equations
   */

  void
  solve_linear_system(bool   initial_step,
                      double relative_residual,
                      double minimum_residual) override;

  /**
   * GMRES solver with ILU preconditioning
   */
  void
  solve_system_GMRES(bool   initial_step,
                     double relative_residual,
                     double minimum_residual);
  /**
   * GMRES solver with AMG preconditioning
   */
  void
  solve_system_AMG(bool   initial_step,
                   double relative_residual,
                   double minimum_residual);


  /**
   * Members
   */
  TrilinosWrappers::BlockSparsityPattern sparsity_pattern;
  TrilinosWrappers::BlockSparseMatrix    system_matrix;
  TrilinosWrappers::SparseMatrix         pressure_mass_matrix;

  std::vector<types::global_dof_index> dofs_per_block;

  const double gamma = 1;
};

template <class BSPreconditioner>
class BlockSchurPreconditioner : public Subscriptor
{
public:
  BlockSchurPreconditioner(double                                     gamma,
                           double                                     viscosity,
                           const TrilinosWrappers::BlockSparseMatrix &S,
                           const TrilinosWrappers::SparseMatrix &     P,
                           const BSPreconditioner * p_amat_preconditioner,
                           const BSPreconditioner * p_pmass_preconditioner,
                           Parameters::LinearSolver solver_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const double                               gamma;
  const double                               viscosity;
  const Parameters::LinearSolver             linear_solver_parameters;
  const TrilinosWrappers::BlockSparseMatrix &stokes_matrix;
  const TrilinosWrappers::SparseMatrix &     pressure_mass_matrix;
  const BSPreconditioner *                   amat_preconditioner;
  const BSPreconditioner *                   pmass_preconditioner;
};

// We can notice that the initialization of the inverse of the matrix at the
// top left corner is completed in the constructor. If so, every application
// of the preconditioner then no longer requires the computation of the
// matrix factors.
template <typename BSPreconditioner>
BlockSchurPreconditioner<BSPreconditioner>::BlockSchurPreconditioner(
  double                                     gamma,
  double                                     viscosity,
  const TrilinosWrappers::BlockSparseMatrix &S,
  const TrilinosWrappers::SparseMatrix &     P,
  const BSPreconditioner *                   p_amat_preconditioner,
  const BSPreconditioner *                   p_pmass_preconditioner,

  Parameters::LinearSolver p_solver_parameters)
  : gamma(gamma)
  , viscosity(viscosity)
  , linear_solver_parameters(p_solver_parameters)
  , stokes_matrix(S)
  , pressure_mass_matrix(P)
  , amat_preconditioner(p_amat_preconditioner)
  , pmass_preconditioner(p_pmass_preconditioner)
{}

template <class BSPreconditioner>
void
BlockSchurPreconditioner<BSPreconditioner>::vmult(
  TrilinosWrappers::MPI::BlockVector &      dst,
  const TrilinosWrappers::MPI::BlockVector &src) const
{
  //    TimerOutput computing_timer(std::cout,
  //                                TimerOutput::summary,
  //                                TimerOutput::wall_times);

  TrilinosWrappers::MPI::Vector utmp(src.block(0));
  {
    //        computing_timer.enter_section("Pressure");
    SolverControl solver_control(
      linear_solver_parameters.max_iterations,
      std::max(
        1e-3 * src.block(1).l2_norm(),
        // linear_solver_parameters.relative_residual*src.block(0).l2_norm(),
        linear_solver_parameters.minimum_residual));
    TrilinosWrappers::SolverCG cg(solver_control);

    dst.block(1) = 0.0;
    cg.solve(pressure_mass_matrix,
             dst.block(1),
             src.block(1),
             *pmass_preconditioner);
    dst.block(1) *= -(viscosity + gamma);
    //        computing_timer.exit_section("Pressure");
  }

  {
    //        computing_timer.enter_section("Operations");
    stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
    utmp *= -1.0;
    utmp += src.block(0);
    //        computing_timer.exit_section("Operations");
  }
  {
    //        computing_timer.enter_section("A Matrix");
    SolverControl solver_control(
      linear_solver_parameters.max_iterations,
      std::max(1e-1 * src.block(0).l2_norm(),
               linear_solver_parameters.minimum_residual));



    // TrilinosWrappers::SolverBicgstab solver(solver_control);
    TrilinosWrappers::SolverGMRES solver(solver_control);

    // A_inverse.solve(stokes_matrix.block(0, 0),dst.block(0), utmp,
    // mp_preconditioner);
    solver.solve(stokes_matrix.block(0, 0),
                 dst.block(0),
                 utmp,
                 *amat_preconditioner);
    //        computing_timer.exit_section("A Matrix");
  }
}



template <class PreconditionerMp>
class BlockDiagPreconditioner : public Subscriptor
{
public:
  BlockDiagPreconditioner(const TrilinosWrappers::BlockSparseMatrix &S,
                          const PreconditionerMp &Mppreconditioner,
                          SolverControl &         A_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const TrilinosWrappers::BlockSparseMatrix & stokes_matrix;
  TrilinosWrappers::PreconditionILU           amat_preconditioner;
  TrilinosWrappers::PreconditionILU           pmat_preconditioner;
  const PreconditionerMp &                    mp_preconditioner;
  SolverFGMRES<TrilinosWrappers::MPI::Vector> A_inverse;
};

template <class PreconditionerMp>
BlockDiagPreconditioner<PreconditionerMp>::BlockDiagPreconditioner(
  const TrilinosWrappers::BlockSparseMatrix &S,
  const PreconditionerMp &                   Mppreconditioner,
  SolverControl &                            A_parameters)
  : stokes_matrix(S)
  , mp_preconditioner(Mppreconditioner)
  , A_inverse(A_parameters)
{
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(0,
                                                                          1e-12,
                                                                          1,
                                                                          0);
  amat_preconditioner.initialize(stokes_matrix.block(0, 0),
                                 preconditionerOptions);

  TrilinosWrappers::PreconditionILU::AdditionalData pmat_preconditionerOptions(
    0, 1e-12, 1, 0);
  TrilinosWrappers::PreconditionILU pmat_preconditioner;
  pmat_preconditioner.initialize(stokes_matrix.block(1, 1),
                                 pmat_preconditionerOptions);
}

template <class PreconditionerMp>
void
BlockDiagPreconditioner<PreconditionerMp>::vmult(
  TrilinosWrappers::MPI::BlockVector &      dst,
  const TrilinosWrappers::MPI::BlockVector &src) const
{
  {
    SolverControl              solver_control(100000, 1e-12);
    TrilinosWrappers::SolverCG cg(solver_control);

    cg.solve(stokes_matrix.block(1, 1),
             dst.block(1),
             src.block(1),
             pmat_preconditioner);
  }
  {
    SolverControl solver_control(10000, 1e-12);

    TrilinosWrappers::SolverGMRES solver(solver_control);
    solver.solve(stokes_matrix.block(0, 0),
                 dst.block(0),
                 src.block(0),
                 amat_preconditioner);
  }
}


// Constructor for class GDNavierStokesSolver
template <int dim>
GDNavierStokesSolver<dim>::GDNavierStokesSolver(
  NavierStokesSolverParameters<dim> &p_nsparam,
  const unsigned int                 degreeVelocity,
  const unsigned int                 degreePressure)
  : NavierStokesBase<dim,
                     TrilinosWrappers::MPI::BlockVector,
                     std::vector<IndexSet>>(p_nsparam,
                                            degreeVelocity,
                                            degreePressure)
{}

template <int dim>
GDNavierStokesSolver<dim>::~GDNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
GDNavierStokesSolver<dim>::set_solution_vector(double value)
{
  this->present_solution = value;
}

template <int dim>
void
GDNavierStokesSolver<dim>::setup_dofs()
{
  TimerOutput::Scope t(this->computing_timer, "setup_dofs");

  system_matrix.clear();

  this->dof_handler.distribute_dofs(this->fe);
  // DoFRenumbering::Cuthill_McKee(this->dof_handler);


  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(this->dof_handler, block_component);
  dofs_per_block.resize(2);
  DoFTools::count_dofs_per_block(this->dof_handler,
                                 dofs_per_block,
                                 block_component);
  unsigned int dof_u = dofs_per_block[0];
  unsigned int dof_p = dofs_per_block[1];

  this->locally_owned_dofs.resize(2);
  this->locally_owned_dofs[0] =
    this->dof_handler.locally_owned_dofs().get_view(0, dof_u);
  this->locally_owned_dofs[1] =
    this->dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

  IndexSet locally_relevant_dofs_acquisition;
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          locally_relevant_dofs_acquisition);
  this->locally_relevant_dofs.resize(2);
  this->locally_relevant_dofs[0] =
    locally_relevant_dofs_acquisition.get_view(0, dof_u);
  this->locally_relevant_dofs[1] =
    locally_relevant_dofs_acquisition.get_view(dof_u, dof_u + dof_p);

  const MappingQ<dim>        mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  {
    this->nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->nonzero_constraints);
    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundaryConditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundaryConditions.type[i_bc] ==
            BoundaryConditions::noslip)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundaryConditions.id[i_bc],
              ZeroFunction<dim>(dim + 1),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }
        else if (this->nsparam.boundaryConditions.type[i_bc] ==
                 BoundaryConditions::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundaryConditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->nonzero_constraints);
          }
        else if (this->nsparam.boundaryConditions.type[i_bc] ==
                 BoundaryConditions::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundaryConditions.id[i_bc],
              FunctionDefined<dim>(
                &this->nsparam.boundaryConditions.bcFunctions[i_bc].u,
                &this->nsparam.boundaryConditions.bcFunctions[i_bc].v,
                &this->nsparam.boundaryConditions.bcFunctions[i_bc].w),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }

        else if (this->nsparam.boundaryConditions.type[i_bc] ==
                 BoundaryConditions::periodic)
          {
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
              this->dof_handler,
              this->nsparam.boundaryConditions.id[i_bc],
              this->nsparam.boundaryConditions.periodic_id[i_bc],
              this->nsparam.boundaryConditions.periodic_direction[i_bc],
              this->nonzero_constraints);
          }
      }
  }
  this->nonzero_constraints.close();

  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);

    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundaryConditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundaryConditions.type[i_bc] ==
            BoundaryConditions::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundaryConditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->zero_constraints);
          }
        else if (this->nsparam.boundaryConditions.type[i_bc] ==
                 BoundaryConditions::periodic)
          {
            DoFTools::make_periodicity_constraints<DoFHandler<dim>>(
              this->dof_handler,
              this->nsparam.boundaryConditions.id[i_bc],
              this->nsparam.boundaryConditions.periodic_id[i_bc],
              this->nsparam.boundaryConditions.periodic_direction[i_bc],
              this->zero_constraints);
          }
        else // if(nsparam.boundaryConditions.boundaries[i_bc].type==Parameters::noslip
             // || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundaryConditions.id[i_bc],
              ZeroFunction<dim>(dim + 1),
              this->zero_constraints,
              this->fe.component_mask(velocities));
          }
      }
  }
  this->zero_constraints.close();

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  this->solution_m1.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m2.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m3.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
                           this->mpi_communicator);

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);


  sparsity_pattern.reinit(this->locally_owned_dofs,
                          this->locally_owned_dofs,
                          this->locally_relevant_dofs,
                          MPI_COMM_WORLD);

  Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
  for (unsigned int c = 0; c < dim + 1; ++c)
    for (unsigned int d = 0; d < dim + 1; ++d)
      if (!((c == dim) && (d == dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::always;

  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  coupling,
                                  sparsity_pattern,
                                  this->nonzero_constraints,
                                  true,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));

  sparsity_pattern.compress();

  system_matrix.reinit(sparsity_pattern);
  pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));


  this->globalVolume_ = GridTools::volume(*this->triangulation);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << this->globalVolume_
              << std::endl;
}

template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme>
void
GDNavierStokesSolver<dim>::assembleGD()
{
  double viscosity = this->nsparam.physicalProperties.viscosity;

  Function<dim> *l_forcing_function = this->forcing_function;

  if (assemble_matrix)
    system_matrix = 0;

  this->system_rhs = 0;

  QGauss<dim>         quadrature_formula(this->degreeQuadrature_);
  const MappingQ<dim> mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);

  FEValues<dim> fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));


  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // For the linearized system, we create temporary storage for present
  // velocity and gradient, and present pressure. In practice, they are all
  // obtained through their shape functions at quadrature points.

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
  std::vector<double>         present_pressure_values(n_q_points);

  std::vector<double>         div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  Tensor<1, dim> force;

  // Get the BDF coefficients
  Vector<double> alpha_bdf;

  if (scheme == Parameters::SimulationControl::bdf1)
    alpha_bdf = bdf_coefficients(1, this->simulationControl.getTimeSteps());

  if (scheme == Parameters::SimulationControl::bdf2)
    alpha_bdf = bdf_coefficients(2, this->simulationControl.getTimeSteps());

  if (scheme == Parameters::SimulationControl::bdf3)
    alpha_bdf = bdf_coefficients(3, this->simulationControl.getTimeSteps());

  // Values at previous time step for backward Euler scheme
  std::vector<Tensor<1, dim>> p1_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p2_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p3_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> p4_velocity_values(n_q_points);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          local_matrix = 0;
          local_rhs    = 0;

          fe_values[velocities].get_function_values(this->evaluation_point,
                                                    present_velocity_values);

          fe_values[velocities].get_function_gradients(
            this->evaluation_point, present_velocity_gradients);

          fe_values[pressure].get_function_values(this->evaluation_point,
                                                  present_pressure_values);

          if (scheme != Parameters::SimulationControl::steady)
            fe_values[velocities].get_function_values(this->solution_m1,
                                                      p1_velocity_values);

          if (scheme == Parameters::SimulationControl::bdf2 ||
              scheme == Parameters::SimulationControl::bdf3)
            fe_values[velocities].get_function_values(this->solution_m2,
                                                      p2_velocity_values);

          if (scheme == Parameters::SimulationControl::bdf3)
            fe_values[velocities].get_function_values(this->solution_m3,
                                                      p3_velocity_values);

          if (l_forcing_function)
            l_forcing_function->vector_value_list(
              fe_values.get_quadrature_points(), rhs_force);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Establish the force vector
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  force[i] = rhs_force[q](component_i);
                }

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  phi_u[k]      = fe_values[velocities].value(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  if (assemble_matrix)
                    {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j)
                        {
                          local_matrix(i, j) +=
                            (viscosity *
                               scalar_product(grad_phi_u[j], grad_phi_u[i]) +
                             present_velocity_gradients[q] * phi_u[j] *
                               phi_u[i] +
                             grad_phi_u[j] * present_velocity_values[q] *
                               phi_u[i] -
                             div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] +
                             gamma * div_phi_u[j] * div_phi_u[i] +
                             phi_p[i] * phi_p[j]) *
                            fe_values.JxW(q);

                          // Mass matrix
                          if (scheme == Parameters::SimulationControl::bdf1 ||
                              scheme == Parameters::SimulationControl::bdf2 ||
                              scheme == Parameters::SimulationControl::bdf3)
                            local_matrix(i, j) += phi_u[j] * phi_u[i] *
                                                  alpha_bdf[0] *
                                                  fe_values.JxW(q);
                        }
                    }

                  double present_velocity_divergence =
                    trace(present_velocity_gradients[q]);
                  local_rhs(i) +=
                    (-viscosity * scalar_product(present_velocity_gradients[q],
                                                 grad_phi_u[i]) -
                     present_velocity_gradients[q] *
                       present_velocity_values[q] * phi_u[i] +
                     present_pressure_values[q] * div_phi_u[i] +
                     present_velocity_divergence * phi_p[i] -
                     gamma * present_velocity_divergence * div_phi_u[i] +
                     force * phi_u[i]) *
                    fe_values.JxW(q);

                  if (scheme == Parameters::SimulationControl::bdf1)
                    local_rhs(i) -=
                      alpha_bdf[0] *
                      (present_velocity_values[q] - p1_velocity_values[q]) *
                      phi_u[i] * fe_values.JxW(q);

                  if (scheme == Parameters::SimulationControl::bdf2)
                    local_rhs(i) -=
                      (alpha_bdf[0] * (present_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[1] * (p1_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[2] * (p2_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);

                  if (scheme == Parameters::SimulationControl::bdf3)
                    local_rhs(i) -=
                      (alpha_bdf[0] * (present_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[1] * (p1_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[2] * (p2_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[3] * (p3_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);

          const AffineConstraints<double> &constraints_used =
            this->zero_constraints;

          if (assemble_matrix)
            {
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          local_dof_indices,
                                                          system_matrix,
                                                          this->system_rhs);
            }
          else
            {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          this->system_rhs);
            }
        }
    }

  if (assemble_matrix)
    {
      system_matrix.compress(VectorOperation::add);

      // Finally we move pressure mass matrix into a separate matrix:
      pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
      pressure_mass_matrix.copy_from(system_matrix.block(1, 1));

      // Note that settings this pressure block to zero is not identical to
      // not assembling anything in this block, because this operation here
      // will (incorrectly) delete diagonal entries that come in from
      // hanging node constraints for pressure DoFs. This means that our
      // whole system matrix will have rows that are completely
      // zero. Luckily, FGMRES handles these rows without any problem.
      system_matrix.block(1, 1) = 0;
    }
  this->system_rhs.compress(VectorOperation::add);
}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GDNavierStokesSolver<dim>::set_initial_condition(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::InitialConditionType::L2projection)
    {
      assemble_L2_projection();
      solve_L2_system(true, 1e-15, 1e-15);
      this->present_solution = this->newton_update;
      this->finish_time_step();
      this->postprocess(true);
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      set_nodal_values();
      this->finish_time_step();
      this->postprocess(true);
    }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      set_nodal_values();
      double viscosity = this->nsparam.physicalProperties.viscosity;
      this->nsparam.physicalProperties.viscosity =
        this->nsparam.initialCondition->viscosity;
      Parameters::SimulationControl::TimeSteppingMethod previousControl =
        this->simulationControl.getMethod();
      this->simulationControl.setMethod(Parameters::SimulationControl::steady);
      PhysicsSolver<TrilinosWrappers::MPI::BlockVector>::
        solve_non_linear_system(Parameters::SimulationControl::steady, false);
      this->simulationControl.setMethod(previousControl);
      this->finish_time_step();
      this->postprocess(true);
      this->simulationControl.setMethod(previousControl);
      this->nsparam.physicalProperties.viscosity = viscosity;
    }
  else
    {
      throw std::runtime_error("GDNS - Initial condition could not be set");
    }
}

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  this->system_rhs = 0;
  QGauss<dim>                 quadrature_formula(this->degreeQuadrature_);
  const MappingQ<dim>         mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  FEValues<dim>               fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int          dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int          n_q_points    = quadrature_formula.size();
  FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              local_rhs(dofs_per_cell);
  std::vector<Vector<double>> initial_velocity(n_q_points,
                                               Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  Tensor<1, dim> rhs_initial_velocity_pressure;
  double         rhs_initial_pressure;

  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->nsparam.initialCondition->uvwp.vector_value_list(
            fe_values.get_quadrature_points(), initial_velocity);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_p[k] = fe_values[pressure].value(k, q);
                  phi_u[k] = fe_values[velocities].value(k, q);
                }

              // Establish the rhs tensor operator
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe.system_to_component_index(i).first;
                  rhs_initial_velocity_pressure[i] =
                    initial_velocity[q](component_i);
                }
              rhs_initial_pressure = initial_velocity[q](dim);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_u[j] * phi_u[i]) * fe_values.JxW(q);
                      local_matrix(i, j) +=
                        (phi_p[j] * phi_p[i]) * fe_values.JxW(q);
                    }
                  local_rhs(i) += (phi_u[i] * rhs_initial_velocity_pressure +
                                   phi_p[i] * rhs_initial_pressure) *
                                  fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          const AffineConstraints<double> &constraints_used =
            this->nonzero_constraints;
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      this->system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GDNavierStokesSolver<dim>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  const MappingQ<dim>              mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initialCondition->uvwp,
                           this->newton_update,
                           this->fe.component_mask(velocities));
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initialCondition->uvwp,
                           this->newton_update,
                           this->fe.component_mask(pressure));
  this->nonzero_constraints.distribute(this->newton_update);
  this->present_solution = this->newton_update;
}

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_matrix_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (time_stepping_method == Parameters::SimulationControl::bdf1)
    assembleGD<true, Parameters::SimulationControl::bdf1>();
  else if (time_stepping_method == Parameters::SimulationControl::bdf2)
    assembleGD<true, Parameters::SimulationControl::bdf2>();
  else if (time_stepping_method == Parameters::SimulationControl::bdf3)
    assembleGD<true, Parameters::SimulationControl::bdf3>();
  else if (time_stepping_method == Parameters::SimulationControl::steady)
    assembleGD<true, Parameters::SimulationControl::steady>();
}

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (time_stepping_method == Parameters::SimulationControl::bdf1)
    assembleGD<false, Parameters::SimulationControl::bdf1>();
  else if (time_stepping_method == Parameters::SimulationControl::bdf2)
    assembleGD<false, Parameters::SimulationControl::bdf2>();
  else if (time_stepping_method == Parameters::SimulationControl::bdf3)
    assembleGD<false, Parameters::SimulationControl::bdf3>();
  else if (time_stepping_method == Parameters::SimulationControl::steady)
    assembleGD<false, Parameters::SimulationControl::steady>();
}

template <int dim>
void
GDNavierStokesSolver<dim>::solve_L2_system(const bool initial_step,
                                           double     absolute_residual,
                                           double     relative_residual)
{
  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linearSolver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linearSolver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
  TrilinosWrappers::PreconditionILU                pmass_preconditioner;

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = this->nsparam.linearSolver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linearSolver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linearSolver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  pmass_preconditioner.initialize(pressure_mass_matrix, preconditionerOptions);


  const BlockDiagPreconditioner<TrilinosWrappers::PreconditionILU>
    preconditioner(system_matrix, pmass_preconditioner, solver_control);

  // preconditioner.initialize(system_matrix, preconditionerOptions);

  solver.solve(system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               preconditioner);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}

template <int dim>
void
GDNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                               double     absolute_residual,
                                               double     relative_residual)
{
  if (this->nsparam.linearSolver.solver == this->nsparam.linearSolver.gmres)
    solve_system_GMRES(initial_step, absolute_residual, relative_residual);
  else if (this->nsparam.linearSolver.solver == this->nsparam.linearSolver.amg)
    solve_system_AMG(initial_step, absolute_residual, relative_residual);
  else
    throw(std::runtime_error("This solver is not allowed"));
}

template <int dim>
void
GDNavierStokesSolver<dim>::solve_system_GMRES(const bool initial_step,
                                              double     absolute_residual,
                                              double     relative_residual)
{
  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linearSolver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = this->nsparam.linearSolver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linearSolver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linearSolver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  TrilinosWrappers::PreconditionILU velocity_preconditioner;
  velocity_preconditioner.initialize(system_matrix.block(0, 0),
                                     preconditionerOptions);

  TrilinosWrappers::PreconditionILU pressure_preconditioner;
  pressure_preconditioner.initialize(system_matrix.block(1, 1),
                                     preconditionerOptions);

  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linearSolver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  const BlockSchurPreconditioner<TrilinosWrappers::PreconditionILU>
    preconditioner(gamma,
                   this->nsparam.physicalProperties.viscosity,
                   system_matrix,
                   pressure_mass_matrix,
                   &velocity_preconditioner,
                   &pressure_preconditioner,
                   this->nsparam.linearSolver);

  solver.solve(system_matrix,
               this->newton_update,
               this->system_rhs,
               preconditioner);
  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(this->newton_update);
}

template <int dim>
void
GDNavierStokesSolver<dim>::solve_system_AMG(const bool initial_step,
                                            double     absolute_residual,
                                            double     relative_residual)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linearSolver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }

  //**********************************************
  // Trillinos Wrapper AMG Preconditioner
  //*********************************************
  TrilinosWrappers::PreconditionAMG velocity_preconditioner;
  TrilinosWrappers::PreconditionAMG pressure_preconditioner;

  // Constant modes for velocity
  std::vector<std::vector<bool>> velocity_constant_modes;
  std::vector<bool>              velocity_components(dim + 1, true);
  velocity_components[dim] = false;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   velocity_components,
                                   velocity_constant_modes);

  // Constant modes for pressure
  std::vector<std::vector<bool>> pressure_constant_modes;
  std::vector<bool>              pressure_components(dim + 1, false);
  pressure_components[dim] = true;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   pressure_components,
                                   pressure_constant_modes);

  this->computing_timer.enter_section("AMG_velocity");
  const bool elliptic_velocity     = false;
  bool       higher_order_elements = false;
  if (this->degreeVelocity_ > 1)
    higher_order_elements = true;
  const unsigned int n_cycles = this->nsparam.linearSolver.amg_n_cycles;
  const bool         w_cycle  = this->nsparam.linearSolver.amg_w_cycles;
  const double       aggregation_threshold =
    this->nsparam.linearSolver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->nsparam.linearSolver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->nsparam.linearSolver.amg_smoother_overlap;
  const bool  output_details = false;
  const char *smoother_type  = "Chebyshev";  //"ILU";
  const char *coarse_type    = "Amesos-KLU"; //"ILU";

  TrilinosWrappers::PreconditionAMG::AdditionalData
    velocity_preconditioner_options(elliptic_velocity,
                                    higher_order_elements,
                                    n_cycles,
                                    w_cycle,
                                    aggregation_threshold,
                                    velocity_constant_modes,
                                    smoother_sweeps,
                                    smoother_overlap,
                                    output_details,
                                    smoother_type,
                                    coarse_type);

  Teuchos::ParameterList              velocity_parameter_ml;
  std::unique_ptr<Epetra_MultiVector> velocity_distributed_constant_modes;
  velocity_preconditioner_options.set_parameters(
    velocity_parameter_ml,
    velocity_distributed_constant_modes,
    system_matrix.block(0, 0));
  velocity_preconditioner.initialize(system_matrix.block(0, 0),
                                     velocity_parameter_ml);
  this->computing_timer.exit_section("AMG_velocity");


  this->computing_timer.enter_section("AMG_pressure");
  const bool elliptic_pressure = true;
  higher_order_elements        = false;
  if (this->degreePressure_ > 1)
    higher_order_elements = true;
  TrilinosWrappers::PreconditionAMG::AdditionalData
                                      pressure_preconditioner_options(elliptic_pressure,
                                    higher_order_elements,
                                    n_cycles,
                                    w_cycle,
                                    aggregation_threshold,
                                    pressure_constant_modes,
                                    smoother_sweeps,
                                    smoother_overlap,
                                    output_details,
                                    smoother_type,
                                    coarse_type);
  Teuchos::ParameterList              pressure_parameter_ml;
  std::unique_ptr<Epetra_MultiVector> pressure_distributed_constant_modes;
  velocity_preconditioner_options.set_parameters(
    pressure_parameter_ml,
    pressure_distributed_constant_modes,
    system_matrix.block(0, 0));
  pressure_preconditioner.initialize(system_matrix.block(1, 1),
                                     pressure_parameter_ml);
  this->computing_timer.exit_section("AMG_pressure");



  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linearSolver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  const BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG>
    preconditioner(gamma,
                   this->nsparam.physicalProperties.viscosity,
                   system_matrix,
                   pressure_mass_matrix,
                   &velocity_preconditioner,
                   &pressure_preconditioner,
                   this->nsparam.linearSolver);

  solver.solve(system_matrix,
               this->newton_update,
               this->system_rhs,
               preconditioner);
  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(this->newton_update);
}

/*
 * Generic CFD Solver application
 * Handles the majority of the cases for the GD-NS solver
 */
template <int dim>
void
GDNavierStokesSolver<dim>::solve()
{
  this->read_mesh();
  this->create_manifolds();

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initialCondition->type,
                              this->nsparam.restartParameters.restart);

  while (this->simulationControl.integrate())
    {
      printTime(this->pcout, this->simulationControl);
      if (!this->simulationControl.firstIter())
        {
          NavierStokesBase<dim,
                           TrilinosWrappers::MPI::BlockVector,
                           std::vector<IndexSet>>::refine_mesh();
        }
      this->iterate(this->simulationControl.firstIter());
      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

#endif
