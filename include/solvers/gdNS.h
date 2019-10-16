/* ---------------------------------------------------------------------
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
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
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
#include <core/parameters.h>
#include <core/pvdhandler.h>
#include <core/simulationcontrol.h>

#include "boundary_conditions.h"
#include "manifolds.h"
#include "navier_stokes_base.h"
#include "navier_stokes_solver_parameters.h"
#include "postprocessors.h"

// Std
#include <fstream>
#include <iostream>

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
  : public NavierStokesBase<dim, TrilinosWrappers::MPI::BlockVector>
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
                        time_stepping_method);

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method);

  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme>
  void
  assembleGD();
  void
  assemble_L2_projection();

  void
  refine_mesh();

  void
  refine_mesh_Kelly();

  void
  refine_mesh_uniform();

  void
  set_nodal_values();

  virtual void
  setup_dofs();

  void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool                             restart = false);
  void
  set_solution_vector(double value);

  virtual void
  solve_non_linear_system(
    const Parameters::SimulationControl::TimeSteppingMethod
               time_stepping_method,
    const bool first_iteration)
  {
    newton_iteration(time_stepping_method, first_iteration);
  }

  void
  newton_iteration(const Parameters::SimulationControl::TimeSteppingMethod,
                   const bool is_initial_step);

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
                      double minimum_residual);


  /**
   * Members
   */
  std::vector<IndexSet> locally_owned_dofs;
  std::vector<IndexSet> locally_relevant_dofs;

  TrilinosWrappers::BlockSparsityPattern sparsity_pattern;
  TrilinosWrappers::BlockSparseMatrix    system_matrix;
  TrilinosWrappers::SparseMatrix         pressure_mass_matrix;

  std::vector<types::global_dof_index> dofs_per_block;

  TrilinosWrappers::MPI::BlockVector system_rhs;
  TrilinosWrappers::MPI::BlockVector evaluation_point;
  TrilinosWrappers::MPI::BlockVector local_evaluation_point;

  const double gamma = 1;
};

template <class PreconditionerMp>
class BlockSchurPreconditioner : public Subscriptor
{
public:
  BlockSchurPreconditioner(double                                     gamma,
                           double                                     viscosity,
                           const TrilinosWrappers::BlockSparseMatrix &S,
                           const TrilinosWrappers::SparseMatrix &     P,
                           const PreconditionerMp &Mppreconditioner,
                           SolverControl &         A_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const double                               gamma;
  const double                               viscosity;
  const TrilinosWrappers::BlockSparseMatrix &stokes_matrix;
  const TrilinosWrappers::SparseMatrix &     pressure_mass_matrix;
  TrilinosWrappers::PreconditionILU          amat_preconditioner;
  TrilinosWrappers::PreconditionILU          pmass_preconditioner;

  const PreconditionerMp &                    mp_preconditioner;
  SolverFGMRES<TrilinosWrappers::MPI::Vector> A_inverse;
};

// We can notice that the initialization of the inverse of the matrix at the
// top left corner is completed in the constructor. If so, every application
// of the preconditioner then no longer requires the computation of the
// matrix factors.
template <class PreconditionerMp>
BlockSchurPreconditioner<PreconditionerMp>::BlockSchurPreconditioner(
  double                                     gamma,
  double                                     viscosity,
  const TrilinosWrappers::BlockSparseMatrix &S,
  const TrilinosWrappers::SparseMatrix &     P,
  const PreconditionerMp &                   Mppreconditioner,
  SolverControl &                            A_parameters)
  : gamma(gamma)
  , viscosity(viscosity)
  , stokes_matrix(S)
  , pressure_mass_matrix(P)
  , mp_preconditioner(Mppreconditioner)
  , A_inverse(A_parameters)
{
  TrilinosWrappers::PreconditionILU::AdditionalData amat_preconditionerOptions(
    0, 1e-10, 1, 0);
  amat_preconditioner.initialize(stokes_matrix.block(0, 0),
                                 amat_preconditionerOptions);

  TrilinosWrappers::PreconditionILU::AdditionalData pmass_preconditionerOptions(
    0, 1e-10, 1, 0);
  TrilinosWrappers::PreconditionILU pmass_preconditioner;
  pmass_preconditioner.initialize(pressure_mass_matrix,
                                  pmass_preconditionerOptions);
}

template <class PreconditionerMp>
void
BlockSchurPreconditioner<PreconditionerMp>::vmult(
  TrilinosWrappers::MPI::BlockVector &      dst,
  const TrilinosWrappers::MPI::BlockVector &src) const
{
  //  MPI_Comm   this->mpi_communicator(MPI_COMM_WORLD);
  //  ConditionalOStream pcout(std::cout,
  //  (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0));
  //  TimerOutput computing_timer(this->mpi_communicator,
  //                              pcout,
  //                              TimerOutput::summary,
  //                              TimerOutput::wall_times);

  TrilinosWrappers::MPI::Vector utmp(src.block(0));
  {
    //    computing_timer.enter_section("Pressure");
    SolverControl              solver_control(100000,
                                 std::max(1e-3 * src.block(0).l2_norm(),
                                          1e-12));
    TrilinosWrappers::SolverCG cg(solver_control);

    dst.block(1) = 0.0;
    cg.solve(pressure_mass_matrix,
             dst.block(1),
             src.block(1),
             pmass_preconditioner);
    dst.block(1) *= -(viscosity + gamma);
    //    computing_timer.exit_section("Pressure");
  }

  {
    //    computing_timer.enter_section("Operations");
    stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
    utmp *= -1.0;
    utmp += src.block(0);
    //    computing_timer.exit_section("Operations");
  }
  {
    //    computing_timer.enter_section("A Matrix");
    SolverControl solver_control(10000,
                                 std::max(1e-2 * src.block(0).l2_norm(),
                                          1e-12));



    TrilinosWrappers::SolverBicgstab solver(solver_control);

    // A_inverse.solve(stokes_matrix.block(0, 0),dst.block(0), utmp,
    // mp_preconditioner);
    solver.solve(stokes_matrix.block(0, 0),
                 dst.block(0),
                 utmp,
                 amat_preconditioner);
    //    computing_timer.exit_section("A Matrix");
  }
}

template <class PreconditionerMp>
class BlockCGPreconditioner : public Subscriptor
{
public:
  BlockCGPreconditioner(double                                     gamma,
                        double                                     viscosity,
                        const TrilinosWrappers::BlockSparseMatrix &S,
                        const PreconditionerMp &Mppreconditioner,
                        SolverControl &         A_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const double                                gamma;
  const double                                viscosity;
  const TrilinosWrappers::BlockSparseMatrix & stokes_matrix;
  TrilinosWrappers::PreconditionILU           amat_preconditioner;
  TrilinosWrappers::PreconditionILU           pmat_preconditioner;
  const PreconditionerMp &                    mp_preconditioner;
  SolverFGMRES<TrilinosWrappers::MPI::Vector> A_inverse;
};

template <class PreconditionerMp>
BlockCGPreconditioner<PreconditionerMp>::BlockCGPreconditioner(
  double                                     gamma,
  double                                     viscosity,
  const TrilinosWrappers::BlockSparseMatrix &S,
  const PreconditionerMp &                   Mppreconditioner,
  SolverControl &                            A_parameters)
  : gamma(gamma)
  , viscosity(viscosity)
  , stokes_matrix(S)
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
BlockCGPreconditioner<PreconditionerMp>::vmult(
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
  : NavierStokesBase<dim, TrilinosWrappers::MPI::BlockVector>(p_nsparam,
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

  locally_owned_dofs.resize(2);
  locally_owned_dofs[0] =
    this->dof_handler.locally_owned_dofs().get_view(0, dof_u);
  locally_owned_dofs[1] =
    this->dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

  IndexSet locally_relevant_dofs_acquisition;
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          locally_relevant_dofs_acquisition);
  locally_relevant_dofs.resize(2);
  locally_relevant_dofs[0] =
    locally_relevant_dofs_acquisition.get_view(0, dof_u);
  locally_relevant_dofs[1] =
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

  this->present_solution.reinit(locally_owned_dofs,
                                locally_relevant_dofs,
                                this->mpi_communicator);

  this->solution_m1.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m2.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           this->mpi_communicator);
  this->solution_m3.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           this->mpi_communicator);

  this->newton_update.reinit(locally_owned_dofs, this->mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, this->mpi_communicator);
  local_evaluation_point.reinit(locally_owned_dofs, this->mpi_communicator);


  sparsity_pattern.reinit(locally_owned_dofs,
                          locally_owned_dofs,
                          locally_relevant_dofs,
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


  this->globalVolume_ = GridTools::volume(this->triangulation);

  this->pcout << "   Number of active cells:       "
              << this->triangulation.n_global_active_cells() << std::endl
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

  system_rhs = 0;

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

          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);

          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);

          fe_values[pressure].get_function_values(evaluation_point,
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
                                                          system_rhs);
            }
          else
            {
              constraints_used.distribute_local_to_global(local_rhs,
                                                          local_dof_indices,
                                                          system_rhs);
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
  system_rhs.compress(VectorOperation::add);
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
      newton_iteration(Parameters::SimulationControl::steady, false);
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
  system_matrix = 0;
  system_rhs    = 0;
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
                                                      system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
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
GDNavierStokesSolver<dim>::refine_mesh()
{
  if (this->simulationControl.getIter() %
        this->nsparam.meshAdaptation.frequency ==
      0)
    {
      if (this->nsparam.meshAdaptation.type ==
          this->nsparam.meshAdaptation.kelly)
        refine_mesh_Kelly();
      if (this->nsparam.meshAdaptation.type ==
          this->nsparam.meshAdaptation.uniform)
        refine_mesh_uniform();
    }
}

template <int dim>
void
GDNavierStokesSolver<dim>::refine_mesh_Kelly()
{
  // Time monitoring
  TimerOutput::Scope t(this->computing_timer, "refine");

  Vector<float> estimated_error_per_cell(this->triangulation.n_active_cells());
  const MappingQ<dim>              mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  if (this->nsparam.meshAdaptation.variable ==
      Parameters::MeshAdaptation::pressure)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->degreeQuadrature_ + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(pressure));
    }
  else if (this->nsparam.meshAdaptation.variable ==
           Parameters::MeshAdaptation::velocity)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->degreeQuadrature_ + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(velocity));
    }

  if (this->nsparam.meshAdaptation.fractionType ==
      Parameters::MeshAdaptation::number)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      this->triangulation,
      estimated_error_per_cell,
      this->nsparam.meshAdaptation.fractionRefinement,
      this->nsparam.meshAdaptation.fractionCoarsening,
      this->nsparam.meshAdaptation.maxNbElements);

  else if (this->nsparam.meshAdaptation.fractionType ==
           Parameters::MeshAdaptation::fraction)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      this->triangulation,
      estimated_error_per_cell,
      this->nsparam.meshAdaptation.fractionRefinement,
      this->nsparam.meshAdaptation.fractionCoarsening);

  if (this->triangulation.n_levels() > this->nsparam.meshAdaptation.maxRefLevel)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           this->triangulation.begin_active(
             this->nsparam.meshAdaptation.maxRefLevel);
         cell != this->triangulation.end();
         ++cell)
      cell->clear_refine_flag();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         this->triangulation.begin_active(
           this->nsparam.meshAdaptation.minRefLevel);
       cell !=
       this->triangulation.end_active(this->nsparam.meshAdaptation.minRefLevel);
       ++cell)
    cell->clear_coarsen_flag();

  this->triangulation.prepare_coarsening_and_refinement();

  // Solution transfer objects for all the solutions
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m1(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m2(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m3(this->dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(this->solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(this->solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(this->solution_m3);

  this->triangulation.execute_coarsening_and_refinement();
  setup_dofs();

  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::BlockVector tmp(locally_owned_dofs,
                                         this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m1(locally_owned_dofs,
                                            this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m2(locally_owned_dofs,
                                            this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m3(locally_owned_dofs,
                                            this->mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);
  solution_transfer_m1.interpolate(tmp_m1);
  solution_transfer_m2.interpolate(tmp_m2);
  solution_transfer_m3.interpolate(tmp_m3);

  // Distribute constraints
  this->nonzero_constraints.distribute(tmp);
  this->nonzero_constraints.distribute(tmp_m1);
  this->nonzero_constraints.distribute(tmp_m2);
  this->nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  this->present_solution = tmp;
  this->solution_m1      = tmp_m1;
  this->solution_m2      = tmp_m2;
  this->solution_m3      = tmp_m3;
}

template <int dim>
void
GDNavierStokesSolver<dim>::refine_mesh_uniform()
{
  TimerOutput::Scope t(this->computing_timer, "refine");

  // Solution transfer objects for all the solutions
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m1(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m2(this->dof_handler);
  parallel::distributed::SolutionTransfer<dim,
                                          TrilinosWrappers::MPI::BlockVector>
    solution_transfer_m3(this->dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(this->solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(this->solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(this->solution_m3);

  // Refine
  this->triangulation.refine_global(1);

  setup_dofs();

  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::BlockVector tmp(locally_owned_dofs,
                                         this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m1(locally_owned_dofs,
                                            this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m2(locally_owned_dofs,
                                            this->mpi_communicator);
  TrilinosWrappers::MPI::BlockVector tmp_m3(locally_owned_dofs,
                                            this->mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);
  solution_transfer_m1.interpolate(tmp_m1);
  solution_transfer_m2.interpolate(tmp_m2);
  solution_transfer_m3.interpolate(tmp_m3);

  // Distribute constraints
  this->nonzero_constraints.distribute(tmp);
  this->nonzero_constraints.distribute(tmp_m1);
  this->nonzero_constraints.distribute(tmp_m2);
  this->nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  this->present_solution = tmp;
  this->solution_m1      = tmp_m1;
  this->solution_m2      = tmp_m2;
  this->solution_m3      = tmp_m3;
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
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linearSolver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    locally_owned_dofs, this->mpi_communicator);

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


  const BlockCGPreconditioner<TrilinosWrappers::PreconditionILU> preconditioner(
    gamma,
    this->nsparam.physicalProperties.viscosity,
    system_matrix,
    pmass_preconditioner,
    solver_control);

  // preconditioner.initialize(system_matrix, preconditionerOptions);

  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
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
  TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << std::setprecision(
                       this->nsparam.linearSolver.residual_precision)
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    locally_owned_dofs, this->mpi_communicator);

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

  const BlockSchurPreconditioner<TrilinosWrappers::PreconditionILU>
    preconditioner(gamma,
                   this->nsparam.physicalProperties.viscosity,
                   system_matrix,
                   pressure_mass_matrix,
                   pmass_preconditioner,
                   solver_control);

  solver.solve(system_matrix, this->newton_update, system_rhs, preconditioner);
  if (this->nsparam.linearSolver.verbosity != Parameters::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(this->newton_update);
}


template <int dim>
void
GDNavierStokesSolver<dim>::newton_iteration(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step)
{
  double current_res;
  double last_res;
  bool   first_step = is_initial_step;
  {
    unsigned int outer_iteration = 0;
    last_res                     = 1.0;
    current_res                  = 1.0;
    while ((current_res > this->nsparam.nonLinearSolver.tolerance) &&
           outer_iteration < this->nsparam.nonLinearSolver.maxIterations)
      {
        evaluation_point = this->present_solution;
        assemble_matrix_rhs(time_stepping_method);
        if (outer_iteration == 0)
          {
            current_res = system_rhs.l2_norm();
            last_res    = current_res;
          }
        if (this->nsparam.nonLinearSolver.verbosity != Parameters::quiet)
          this->pcout << "Newton iteration: " << outer_iteration
                      << "  - Residual:  " << current_res << std::endl;
        solve_linear_system(first_step,
                            this->nsparam.linearSolver.minimum_residual,
                            this->nsparam.linearSolver.relative_residual);

        for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
          {
            local_evaluation_point = this->present_solution;
            local_evaluation_point.add(alpha, this->newton_update);
            this->nonzero_constraints.distribute(local_evaluation_point);
            evaluation_point = local_evaluation_point;
            assemble_rhs(time_stepping_method);
            current_res = system_rhs.l2_norm();
            if (this->nsparam.nonLinearSolver.verbosity != Parameters::quiet)
              this->pcout << "\t\talpha = " << std::setw(6) << alpha
                          << std::setw(0) << " res = "
                          << std::setprecision(
                               this->nsparam.nonLinearSolver.display_precision)
                          << current_res << std::endl;
            if (current_res < 0.9 * last_res ||
                last_res < this->nsparam.nonLinearSolver.tolerance)
              break;
          }
        this->present_solution = evaluation_point;
        last_res               = current_res;
        ++outer_iteration;
      }
  }
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
          this->refine_mesh();
        }
      this->iterate(this->simulationControl.firstIter());
      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

#endif
