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

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/gd_navier_stokes.h>
#include <solvers/isothermal_compressible_navier_stokes_vof_assembler.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_vof_assemblers.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/vector_tools.h>



// Constructor for class GDNavierStokesSolver
template <int dim>
GDNavierStokesSolver<dim>::GDNavierStokesSolver(
  SimulationParameters<dim> &p_nsparam)
  : NavierStokesBase<dim,
                     TrilinosWrappers::MPI::BlockVector,
                     std::vector<IndexSet>>(p_nsparam)
{}

template <int dim>
GDNavierStokesSolver<dim>::~GDNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
GDNavierStokesSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Buoyant force
  if (this->simulation_parameters.multiphysics.buoyancy_force)
    {
      this->assemblers.push_back(
        std::make_shared<BuoyancyAssembly<dim>>(this->simulation_control));
    }

  if (this->simulation_parameters.multiphysics.VOF)
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()) &&
          this->simulation_parameters.physical_properties_manager
            .density_is_constant())
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }
      else if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<
              GLSIsothermalCompressibleNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }
      if (!this->simulation_parameters.physical_properties_manager
             .density_is_constant())
        {
          this->assemblers.push_back(
            std::make_shared<
              GLSIsothermalCompressibleNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
      else
        {
          // Core assembler
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
    }
  else
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerBDF<dim>>(
              this->simulation_control));
        }

      // Velocity sources term
      if (this->simulation_parameters.velocity_sources.type ==
          Parameters::VelocitySource::VelocitySourceType::srf)
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSRF<dim>>(
              this->simulation_parameters.velocity_sources));
        }

      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GDNavierStokesAssemblerNonNewtonianCore<dim>>(
              this->simulation_control, gamma));
        }
      else
        {
          // Core default assembler
          if ((this->simulation_parameters.stabilization
                 .use_default_stabilization == true) ||
              this->simulation_parameters.stabilization.stabilization ==
                Parameters::Stabilization::NavierStokesStabilization::grad_div)
            this->assemblers.push_back(
              std::make_shared<GDNavierStokesAssemblerCore<dim>>(
                this->simulation_control, gamma));

          else
            throw std::runtime_error(
              "Using the GD solver with a stabilization other than the grad_div "
              "stabilization will lead to an unstable block solver that is unable to converge");
        }
    }
}

template <int dim>
void
GDNavierStokesSolver<dim>::update_multiphysics_time_average_solution()
{
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->multiphysics->set_block_time_average_solution(
        PhysicsID::fluid_dynamics,
        &this->average_velocities->get_average_velocities());
    }
}

template <int dim>
void
GDNavierStokesSolver<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
  // this->simulation_control->set_assembly_method(this->time_stepping_method);

  this->system_matrix = 0;
  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);
    }


  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &GDNavierStokesSolver::assemble_local_system_matrix,
    &GDNavierStokesSolver::copy_local_matrix_to_global_matrix,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));


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



template <int dim>
void
GDNavierStokesSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);
  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();


  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GDNavierStokesSolver<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              system_matrix);
}


template <int dim>
void
GDNavierStokesSolver<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  // this->simulation_control->set_assembly_method(this->time_stepping_method);

  this->system_rhs = 0;
  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }


  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &GDNavierStokesSolver::assemble_local_system_rhs,
    &GDNavierStokesSolver::copy_local_rhs_to_global_rhs,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}


template <int dim>
void
GDNavierStokesSolver<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      typename DoFHandler<dim>::active_cell_iterator temperature_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_ht);

      scratch_data.reinit_heat_transfer(temperature_cell,
                                        *this->multiphysics->get_solution(
                                          PhysicsID::heat_transfer));
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();
  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}



template <int dim>
void
GDNavierStokesSolver<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}



template <int dim>
void
GDNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  auto &system_rhs = this->system_rhs;
  system_rhs       = 0;
  FEValues<dim>               fe_values(*this->mapping,
                          *this->fe,
                          *this->cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int          dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int          n_q_points    = this->cell_quadrature->size();
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

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->simulation_parameters.initial_condition->uvwp.vector_value_list(
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
                    this->fe->system_to_component_index(i).first;
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
          this->nonzero_constraints.distribute_local_to_global(
            local_matrix,
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
GDNavierStokesSolver<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "Setup DOFs");

  system_matrix.clear();

  this->dof_handler.distribute_dofs(*this->fe);
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

  FEValuesExtractors::Vector velocities(0);

  // Non-zero constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  {
    nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->simulation_parameters.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              nonzero_constraints);
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .u,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .v,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .w),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }

        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
              this->simulation_parameters.boundary_conditions
                .periodic_direction[i_bc],
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);

    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->simulation_parameters.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              this->zero_constraints);
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
              this->simulation_parameters.boundary_conditions
                .periodic_direction[i_bc],
              this->zero_constraints);
          }
        else // if(nsparam.boundaryConditions.boundaries[i_bc].type==Parameters::noslip
          // || Parameters::function)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              this->zero_constraints,
              this->fe->component_mask(velocities));
          }
      }
  }
  this->zero_constraints.close();

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  this->evaluation_point.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  // Initialize vector of previous solutions
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      this->mpi_communicator);
    }

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
                                  nonzero_constraints,
                                  true,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));

  sparsity_pattern.compress();

  system_matrix.reinit(sparsity_pattern);
  pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      AssertThrow(this->simulation_parameters.mesh_adaptation.type ==
                    Parameters::MeshAdaptation::Type::none,
                  ExcMessage(
                    "Time-averaging velocities and calculating reynolds "
                    "stresses are currently unavailable for mesh "
                    "adaptation."));

      this->average_velocities->initialize_vectors(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->fe->n_dofs_per_vertex(),
        this->mpi_communicator);

      if (this->simulation_parameters.restart_parameters.checkpoint)
        {
          this->average_velocities->initialize_checkpoint_vectors(
            this->locally_owned_dofs,
            this->locally_relevant_dofs,
            this->mpi_communicator);
        }
    }


  double global_volume =
    GridTools::volume(*this->triangulation, *this->mapping);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;

  this->multiphysics->set_dof_handler(PhysicsID::fluid_dynamics,
                                      &this->dof_handler);
  this->multiphysics->set_block_solution(PhysicsID::fluid_dynamics,
                                         &this->present_solution);
  this->multiphysics->set_block_previous_solutions(PhysicsID::fluid_dynamics,
                                                   &this->previous_solutions);
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
GDNavierStokesSolver<dim>::set_initial_condition_fd(
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
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->finish_time_step();
    }

  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      this->set_nodal_values();
      std::shared_ptr<RheologicalModel> original_viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Temporarily set the rheology to be newtonian with predefined viscosity
      std::shared_ptr<Newtonian> temporary_rheology =
        std::make_shared<Newtonian>(
          this->simulation_parameters.initial_condition->kinematic_viscosity);

      this->simulation_parameters.physical_properties_manager.set_rheology(
        temporary_rheology);


      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::steady);
      PhysicsSolver<
        TrilinosWrappers::MPI::BlockVector>::solve_non_linear_system(false);
      this->finish_time_step();

      this->simulation_parameters.physical_properties_manager.set_rheology(
        original_viscosity_model);
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
  const double absolute_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .relative_residual;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .solver == Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step,
                       absolute_residual,
                       relative_residual,
                       renewed_matrix);
  else
    AssertThrow(this->simulation_parameters.linear_solver
                    .at(PhysicsID::fluid_dynamics)
                    .solver == Parameters::LinearSolver::SolverType::gmres,
                ExcMessage("This linear solver is not allowed."));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
GDNavierStokesSolver<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;

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
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_ilu_preconditioner),
    &(*pressure_ilu_preconditioner),
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics));
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
  FEValuesExtractors::Vector     velocities(0);

  ComponentMask velocity_components = this->fe->component_mask(velocities);
  // std::vector<bool> velocity_components(dim + 1, true);
  // velocity_components[dim] = false;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   velocity_components,
                                   velocity_constant_modes);

  // Constant modes for pressure
  std::vector<std::vector<bool>> pressure_constant_modes;
  FEValuesExtractors::Scalar     pressure(dim);
  ComponentMask pressure_components = this->fe->component_mask(pressure);
  // std::vector<bool> pressure_components(dim + 1, false);
  // pressure_components[dim] = true;
  DoFTools::extract_constant_modes(this->dof_handler,
                                   pressure_components,
                                   pressure_constant_modes);

  this->computing_timer.enter_subsection("AMG_velocity");
  const bool elliptic_velocity     = false;
  bool       higher_order_elements = false;
  if (this->fe->degree > 1)
    higher_order_elements = true;
  const unsigned int n_cycles =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_n_cycles;
  const bool w_cycle =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_w_cycles;
  const double aggregation_threshold =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_smoother_overlap;
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
  Teuchos::ParameterList pressure_parameter_ml;
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
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_amg_preconditioner),
    &(*pressure_amg_preconditioner),
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics));
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

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }


  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::ilu)
    {
      if (renewed_matrix || velocity_ilu_preconditioner == 0 ||
          pressure_ilu_preconditioner == 0 || system_ilu_preconditioner == 0)
        setup_ILU();
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::amg)
    {
      if (renewed_matrix || velocity_amg_preconditioner == 0 ||
          pressure_amg_preconditioner == 0 || system_amg_preconditioner == 0)
        setup_AMG();
    }

  {
    TimerOutput::Scope t(this->computing_timer, "Solve linear system");
    if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .preconditioner == Parameters::LinearSolver::PreconditionerType::ilu)
      solver.solve(this->system_matrix,
                   this->newton_update,
                   this->system_rhs,
                   *system_ilu_preconditioner);
    else if (this->simulation_parameters.linear_solver
               .at(PhysicsID::fluid_dynamics)
               .preconditioner ==
             Parameters::LinearSolver::PreconditionerType::amg)
      solver.solve(this->system_matrix,
                   this->newton_update,
                   this->system_rhs,
                   *system_amg_preconditioner);
    else
      AssertThrow(
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::ilu ||
          this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::amg,
        ExcMessage("This linear solver does not support this preconditioner."));

    if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .verbosity != Parameters::Verbosity::quiet)
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
  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::BlockVector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);
  SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);
  TrilinosWrappers::PreconditionILU                pmass_preconditioner;

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  pmass_preconditioner.initialize(pressure_mass_matrix, preconditionerOptions);


  const BlockDiagPreconditioner<TrilinosWrappers::PreconditionILU>
    preconditioner(system_matrix, pmass_preconditioner, solver_control);

  // preconditioner.initialize(system_matrix, preconditionerOptions);

  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               preconditioner);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
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
  // This is enforced to 1 right now because it does not provide
  // better speed-up than using MPI. This could be eventually changed...
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();
  this->box_refine_mesh();
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);
  this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          for (unsigned int i = 0;
               i <
               this->simulation_parameters.mesh_adaptation.initial_refinement;
               i++)
            NavierStokesBase<dim,
                             TrilinosWrappers::MPI::BlockVector,
                             std::vector<IndexSet>>::refine_mesh();

          this->iterate();
        }
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
