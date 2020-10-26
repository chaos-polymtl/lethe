/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include "solvers/gd_navier_stokes.h"

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/sdirk.h"
#include "core/time_integration_utilities.h"
#include "core/utilities.h"

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
GDNavierStokesSolver<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_system");

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    assembleGD<true, Parameters::SimulationControl::TimeSteppingMethod::bdf1>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    assembleGD<true, Parameters::SimulationControl::TimeSteppingMethod::bdf2>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    assembleGD<true, Parameters::SimulationControl::TimeSteppingMethod::bdf3>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::steady)
    assembleGD<true,
               Parameters::SimulationControl::TimeSteppingMethod::steady>();
}

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  TimerOutput::Scope t(this->computing_timer, "assemble_rhs");

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    assembleGD<false,
               Parameters::SimulationControl::TimeSteppingMethod::bdf1>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    assembleGD<false,
               Parameters::SimulationControl::TimeSteppingMethod::bdf2>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    assembleGD<false,
               Parameters::SimulationControl::TimeSteppingMethod::bdf3>();
  else if (time_stepping_method ==
           Parameters::SimulationControl::TimeSteppingMethod::steady)
    assembleGD<false,
               Parameters::SimulationControl::TimeSteppingMethod::steady>();
}

template <int dim>
template <bool                                              assemble_matrix,
          Parameters::SimulationControl::TimeSteppingMethod scheme>
void
GDNavierStokesSolver<dim>::assembleGD()
{
  double viscosity = this->nsparam.physical_properties.viscosity;

  Function<dim> *l_forcing_function = this->forcing_function;

  if (assemble_matrix)
    system_matrix = 0;

  this->system_rhs = 0;

  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);

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

  std::vector<double> time_steps =
    this->simulation_control->get_time_steps_vector();

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf1)
    alpha_bdf = bdf_coefficients(1, time_steps);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    alpha_bdf = bdf_coefficients(2, time_steps);

  if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    alpha_bdf = bdf_coefficients(3, time_steps);

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

          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            fe_values[velocities].get_function_values(this->solution_m1,
                                                      p1_velocity_values);

          if (scheme ==
                Parameters::SimulationControl::TimeSteppingMethod::bdf2 ||
              scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
            fe_values[velocities].get_function_values(this->solution_m2,
                                                      p2_velocity_values);

          if (scheme == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
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
                          if (scheme == Parameters::SimulationControl::
                                          TimeSteppingMethod::bdf1 ||
                              scheme == Parameters::SimulationControl::
                                          TimeSteppingMethod::bdf2 ||
                              scheme == Parameters::SimulationControl::
                                          TimeSteppingMethod::bdf3)
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

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf1)
                    local_rhs(i) -=
                      alpha_bdf[0] *
                      (present_velocity_values[q] - p1_velocity_values[q]) *
                      phi_u[i] * fe_values.JxW(q);

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    local_rhs(i) -=
                      (alpha_bdf[0] * (present_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[1] * (p1_velocity_values[q] * phi_u[i]) +
                       alpha_bdf[2] * (p2_velocity_values[q] * phi_u[i])) *
                      fe_values.JxW(q);

                  if (scheme ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
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

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  this->system_rhs = 0;
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int  dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int  n_q_points    = quadrature_formula.size();
  FullMatrix<double>  local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs(dofs_per_cell);
  std::vector<Vector<double>>          initial_velocity(n_q_points,
                                                        Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  Tensor<1, dim> rhs_initial_velocity_pressure;
  double         rhs_initial_pressure;

  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->nsparam.initial_condition->uvwp.vector_value_list(
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
GDNavierStokesSolver<dim>::setup_dofs()
{
  TimerOutput::Scope t(this->computing_timer, "setup_dofs");

  system_matrix.clear();

  this->dof_handler.distribute_dofs(this->fe);
  // DoFRenumbering::Cuthill_McKee(this->dof_handler);


  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(this->dof_handler, block_component);


  // To be used to replace the above part once 9.2 release is out and TRAVIS-CI
  // version is updated
  dofs_per_block =
    DoFTools::count_dofs_per_fe_block(this->dof_handler, block_component);

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

  const MappingQ<dim>        mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  {
    this->nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->nonzero_constraints);
    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->nonzero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].u,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].v,
                &this->nsparam.boundary_conditions.bcFunctions[i_bc].w),
              this->nonzero_constraints,
              this->fe.component_mask(velocities));
          }

        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->nonzero_constraints);
          }
      }
  }
  this->nonzero_constraints.close();

  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);

    for (unsigned int i_bc = 0; i_bc < this->nsparam.boundary_conditions.size;
         ++i_bc)
      {
        if (this->nsparam.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->nsparam.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->zero_constraints);
          }
        else if (this->nsparam.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              this->nsparam.boundary_conditions.periodic_id[i_bc],
              this->nsparam.boundary_conditions.periodic_direction[i_bc],
              this->zero_constraints);
          }
        else // if(nsparam.boundaryConditions.boundaries[i_bc].type==Parameters::noslip
          // || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(
              mapping,
              this->dof_handler,
              this->nsparam.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
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


  double global_volume = GridTools::volume(*this->triangulation);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;
}


template <int dim>
void
GDNavierStokesSolver<dim>::set_solution_vector(double value)
{
  this->present_solution = value;
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
      this->set_nodal_values();
      this->finish_time_step();
      this->postprocess(true);
    }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      this->set_nodal_values();
      double viscosity = this->nsparam.physical_properties.viscosity;
      this->nsparam.physical_properties.viscosity =
        this->nsparam.initial_condition->viscosity;
      PhysicsSolver<TrilinosWrappers::MPI::BlockVector>::
        solve_non_linear_system(
          Parameters::SimulationControl::TimeSteppingMethod::steady,
          false,
          true);
      this->finish_time_step();
      this->postprocess(true);
      this->nsparam.physical_properties.viscosity = viscosity;
    }
  else
    {
      throw std::runtime_error("GDNS - Initial condition could not be set");
    }
}

template <int dim>
void
GDNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                               const bool renewed_matrix)
{
  const double absolute_residual = this->nsparam.linear_solver.minimum_residual;
  const double relative_residual =
    this->nsparam.linear_solver.relative_residual;

  if (this->nsparam.linear_solver.solver ==
      Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step,
                       absolute_residual,
                       relative_residual,
                       renewed_matrix);
  else if (this->nsparam.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::amg)
    solve_system_AMG(initial_step,
                     absolute_residual,
                     relative_residual,
                     renewed_matrix);
  else
    throw(std::runtime_error("This solver is not allowed"));
}

template <int dim>
void
GDNavierStokesSolver<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = this->nsparam.linear_solver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linear_solver.ilu_precond_rtol;

  velocity_ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();
  pressure_ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  velocity_ilu_preconditioner->initialize(system_matrix.block(0, 0),
                                          preconditionerOptions);

  pressure_ilu_preconditioner->initialize(system_matrix.block(1, 1),
                                          preconditionerOptions);
  system_ilu_preconditioner = std::make_shared<
    BlockSchurPreconditioner<TrilinosWrappers::PreconditionILU>>(
    gamma,
    this->nsparam.physical_properties.viscosity,
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_ilu_preconditioner),
    &(*pressure_ilu_preconditioner),
    this->nsparam.linear_solver);
}

template <int dim>
void
GDNavierStokesSolver<dim>::setup_AMG()
{
  TimerOutput::Scope t(this->computing_timer, "setup_AMG");

  //**********************************************
  // Trillinos Wrapper AMG Preconditioner
  //*********************************************

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

  this->computing_timer.enter_subsection("AMG_velocity");
  const bool elliptic_velocity     = false;
  bool       higher_order_elements = false;
  if (this->fe.degree > 1)
    higher_order_elements = true;
  const unsigned int n_cycles = this->nsparam.linear_solver.amg_n_cycles;
  const bool         w_cycle  = this->nsparam.linear_solver.amg_w_cycles;
  const double       aggregation_threshold =
    this->nsparam.linear_solver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->nsparam.linear_solver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->nsparam.linear_solver.amg_smoother_overlap;
  const bool  output_details = false;
  const char *smoother_type  = "Chebyshev";  //"ILU";
  const char *coarse_type    = "Amesos-KLU"; //"ILU";

  velocity_amg_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionAMG>();
  pressure_amg_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionAMG>();

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
  velocity_amg_preconditioner->initialize(system_matrix.block(0, 0),
                                          velocity_parameter_ml);
  this->computing_timer.leave_subsection("AMG_velocity");

  this->computing_timer.enter_subsection("AMG_pressure");
  const bool elliptic_pressure = true;
  higher_order_elements        = false;
  if (this->pressure_fem_degree > 1)
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
  pressure_amg_preconditioner->initialize(system_matrix.block(1, 1),
                                          pressure_parameter_ml);
  this->computing_timer.leave_subsection("AMG_pressure");


  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);


  system_amg_preconditioner = std::make_shared<
    BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG>>(
    gamma,
    this->nsparam.physical_properties.viscosity,
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_amg_preconditioner),
    &(*pressure_amg_preconditioner),
    this->nsparam.linear_solver);
}



template <int dim>
void
GDNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                              const double absolute_residual,
                                              const double relative_residual,
                                              const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }


  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  if (renewed_matrix || velocity_ilu_preconditioner == 0 ||
      pressure_ilu_preconditioner == 0 || system_ilu_preconditioner == 0)
    setup_ILU();

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");
    solver.solve(system_matrix,
                 this->newton_update,
                 this->system_rhs,
                 *system_ilu_preconditioner);
    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    constraints_used.distribute(this->newton_update);
  }
}



template <int dim>
void
GDNavierStokesSolver<dim>::solve_system_AMG(const bool   initial_step,
                                            const double absolute_residual,
                                            const double relative_residual,
                                            const bool   renewed_matrix)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  if (renewed_matrix || velocity_amg_preconditioner == 0 ||
      pressure_amg_preconditioner == 0 || system_amg_preconditioner == 0)
    setup_AMG();


  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  {
    TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

    solver.solve(system_matrix,
                 this->newton_update,
                 this->system_rhs,
                 *system_amg_preconditioner);
    if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    constraints_used.distribute(this->newton_update);
  }
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

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->nsparam.linear_solver.max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
  TrilinosWrappers::PreconditionILU                pmass_preconditioner;

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = this->nsparam.linear_solver.ilu_precond_fill;
  const double ilu_atol = this->nsparam.linear_solver.ilu_precond_atol;
  const double ilu_rtol = this->nsparam.linear_solver.ilu_precond_rtol;
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

  if (this->nsparam.linear_solver.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}


/*
 * Generic CFD Solver application
 * Handles the majority of the cases for the GD-NS solver
 */
template <int dim>
void
GDNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(this->triangulation,
                          this->nsparam.mesh,
                          this->nsparam.manifolds_parameters,
                          this->nsparam.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(this->nsparam.initial_condition->type,
                              this->nsparam.restart_parameters.restart);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);
      if (this->simulation_control->is_at_start())
        this->first_iteration();
      else
        {
          NavierStokesBase<dim,
                           TrilinosWrappers::MPI::BlockVector,
                           std::vector<IndexSet>>::refine_mesh();
          this->iterate();
        }
      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GDNavierStokesSolver<2>;
template class GDNavierStokesSolver<3>;
