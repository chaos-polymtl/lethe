#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/vof.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/vector_tools.h>

template <int dim>
void
VolumeOfFluid<dim>::assemble_matrix_and_rhs()
{
  assemble_system_matrix();
  assemble_system_rhs();
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_rhs()
{
  assemble_system_rhs();
}

template <int dim>
void
VolumeOfFluid<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      if (is_sdirk(this->simulation_control->get_assembly_method()))
        throw std::invalid_argument(
          "SDIRK time-stepping scheme is not supported in the VOF solver ");
      else
        this->assemblers.push_back(
          std::make_shared<VOFAssemblerBDF<dim>>(this->simulation_control));
    }

  // Core assembler
  this->assemblers.push_back(std::make_shared<VOFAssemblerCore<dim>>(
    this->simulation_control, this->simulation_parameters.fem_parameters));
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_system_matrix()
{
  this->system_matrix = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = VOFScratchData<dim>(*this->fe,
                                          *this->cell_quadrature,
                                          *this->fs_mapping,
                                          dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VolumeOfFluid::assemble_local_system_matrix,
                  &VolumeOfFluid::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFScratchData<dim> &                                 scratch_data,
  StabilizedMethodsCopyData &                           copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;


  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics));
    }
  else
    {
      scratch_data.reinit_velocity(
        velocity_cell, *multiphysics->get_solution(PhysicsID::fluid_dynamics));
    }
  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
VolumeOfFluid<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsCopyData &copy_data)
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
VolumeOfFluid<dim>::assemble_system_rhs()
{
  // TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = VOFScratchData<dim>(*this->fe,
                                          *this->cell_quadrature,
                                          *this->fs_mapping,
                                          dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VolumeOfFluid::assemble_local_system_rhs,
                  &VolumeOfFluid::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFScratchData<dim> &                                 scratch_data,
  StabilizedMethodsCopyData &                           copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics));
    }
  else
    {
      scratch_data.reinit_velocity(
        velocity_cell, *multiphysics->get_solution(PhysicsID::fluid_dynamics));
    }

  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
VolumeOfFluid<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              system_rhs);
}



template <int dim>
void
VolumeOfFluid<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  data_out.add_data_vector(dof_handler, present_solution, "phase");
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values_fs(*this->fs_mapping,
                             *fe,
                             *this->error_quadrature,
                             update_values | update_gradients |
                               update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int  n_q_points = this->error_quadrature->size();
  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->phase;
  exact_solution.set_time(simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_fs.reinit(cell);
          fe_values_fs.get_function_values(present_solution, q_scalar_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values_fs.get_quadrature_points(),
                                    q_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = q_scalar_values[q];
              double exact = q_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values_fs.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_volume(int fluid_index)
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values_fs(*this->fs_mapping,
                             *fe,
                             *this->error_quadrature,
                             update_values | update_gradients |
                               update_quadrature_points | update_JxW_values);

  const unsigned int  n_q_points = this->error_quadrature->size();
  std::vector<double> q_scalar_values(n_q_points);

  double volume = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_fs.reinit(cell);
          fe_values_fs.get_function_values(present_solution, q_scalar_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (fluid_index)
                {
                  case 0:
                    {
                      if (q_scalar_values[q] < 0.5)
                        volume += fe_values_fs.JxW(q) * q_scalar_values[q];
                      break;
                    }
                  case 1:
                    {
                      if (q_scalar_values[q] > 0.5)
                        volume += fe_values_fs.JxW(q) * q_scalar_values[q];
                      break;
                    }
                }
            }
        }
    }
  volume = Utilities::MPI::sum(volume, mpi_communicator);
  return volume;
}

template <int dim>
void
VolumeOfFluid<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        error_table.omit_column_from_convergence_rate_evaluation("cells");
      else
        error_table.omit_column_from_convergence_rate_evaluation("time");

      error_table.set_scientific("error_phase", true);
      error_table.set_precision("error_phase",
                                simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::finish_time_step()
{
  percolate_time_vectors();
}

template <int dim>
void
VolumeOfFluid<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() &&
      !first_iteration)
    {
      double phase_error = calculate_L2_error();

      if (simulation_control->is_steady())
        {
          error_table.add_value("cells",
                                this->triangulation->n_global_active_cells());
        }
      else
        {
          error_table.add_value("time", simulation_control->get_current_time());
        }
      error_table.add_value("error_phase", phase_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error phase : " << phase_error << std::endl;
        }
    }

  if (simulation_parameters.multiphysics.conservation_monitoring)
    {
      double volume =
        calculate_volume(simulation_parameters.multiphysics.fluid_index);

      auto         mpi_communicator = triangulation->get_communicator();
      unsigned int this_mpi_process(
        Utilities::MPI::this_mpi_process(mpi_communicator));

      if (this_mpi_process == 0)
        {
          // Set conservation monitoring table
          if (simulation_control->is_steady())
            {
              volume_table_fs.add_value(
                "cells", this->triangulation->n_global_active_cells());
            }
          else
            {
              volume_table_fs.add_value("time",
                                        simulation_control->get_current_time());
            }
          volume_table_fs.add_value("fluid_volume", volume);
          volume_table_fs.set_scientific("fluid_volume", true);

          // Save table to free_surface_monitoring.dat
          std::string   filename = "free_surface_monitoring.dat";
          std::ofstream output(filename.c_str());
          volume_table_fs.write_text(output);
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::modify_solution()
{
  if (this->simulation_parameters.multiphysics.interface_sharpening)
    {
      // Limit the phase fractions between 0 and 1
      update_solution_and_constraints(present_solution);
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        update_solution_and_constraints(previous_solutions[p]);

      // Interface sharpening is done at a constant frequency
      if (this->simulation_control->get_step_number() %
            this->simulation_parameters.interface_sharpening
              .sharpening_frequency ==
          0)
        {
          if (simulation_parameters.linear_solver.verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "Sharpening interface at step "
                          << this->simulation_control->get_step_number()
                          << std::endl;
            }

          // Sharpen the interface of all solutions:
          {
            // Assemble matrix and solve the system for interface sharpening
            assemble_L2_projection_interface_sharpening(present_solution);
            solve_interface_sharpening(present_solution);

            for (unsigned int p = 0; p < previous_solutions.size(); ++p)
              {
                assemble_L2_projection_interface_sharpening(
                  previous_solutions[p]);
                solve_interface_sharpening(previous_solutions[p]);
              }
          }

          // Re limit the phase fractions between 0 and 1 after interface
          // sharpening
          update_solution_and_constraints(present_solution);
          for (unsigned int p = 0; p < previous_solutions.size(); ++p)
            update_solution_and_constraints(previous_solutions[p]);
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::pre_mesh_adaptation()
{
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_communicator();

  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);

  // Distribute constraints
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  // Transfer previous solutions
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      TrilinosWrappers::MPI::Vector tmp_previous_solution(locally_owned_dofs,
                                                          mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }
}

template <int dim>
void
VolumeOfFluid<dim>::compute_kelly(
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (this->simulation_parameters.mesh_adaptation.variable ==
      Parameters::MeshAdaptation::Variable::phase)
    {
      const FEValuesExtractors::Scalar phase(0);

      KellyErrorEstimator<dim>::estimate(
        *this->fs_mapping,
        dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        fe->component_mask(phase));
    }
}

template <int dim>
void
VolumeOfFluid<dim>::write_checkpoint()
{
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;

  sol_set_transfer.push_back(&present_solution);
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&previous_solutions[i]);
    }
  solution_transfer.prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
VolumeOfFluid<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading VOF checkpoint" << std::endl;

  std::vector<TrilinosWrappers::MPI::Vector *> input_vectors(
    1 + previous_solutions.size());
  TrilinosWrappers::MPI::Vector distributed_system(locally_owned_dofs,
                                                   mpi_communicator);
  input_vectors[0] = &distributed_system;


  std::vector<TrilinosWrappers::MPI::Vector> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(locally_owned_dofs, mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  solution_transfer.deserialize(input_vectors);

  present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }
}


template <int dim>
void
VolumeOfFluid<dim>::setup_dofs()
{
  dof_handler.distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  auto mpi_communicator = triangulation->get_communicator();

  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  present_solution.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);

  // Previous solutions for transient schemes
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
    }

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  newton_update.reinit(locally_owned_dofs, mpi_communicator);

  local_evaluation_point.reinit(this->locally_owned_dofs, mpi_communicator);

  {
    nonzero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
  }
  nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            zero_constraints);
  }
  zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(this->dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  system_matrix_phase_fraction.reinit(locally_owned_dofs,
                                      locally_owned_dofs,
                                      dsp,
                                      mpi_communicator);

  complete_system_matrix_phase_fraction.reinit(locally_owned_dofs,
                                               locally_owned_dofs,
                                               dsp,
                                               mpi_communicator);

  nodal_phase_fraction_owned.reinit(locally_owned_dofs, mpi_communicator);

  system_rhs_phase_fraction.reinit(locally_owned_dofs, mpi_communicator);

  complete_system_rhs_phase_fraction.reinit(locally_owned_dofs,
                                            mpi_communicator);

  active_set.clear();
  active_set.set_size(dof_handler.n_dofs());

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of VOF degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

  // Provide the VOF dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::VOF, &dof_handler);
  multiphysics->set_solution(PhysicsID::VOF, &present_solution);

  // the fluid at present iteration is solved BEFORE the VOF (see map
  // solve_pre_fluid defined in multiphysics_interface.h), and after percolate
  // is called for the previous iteration.
  // NB: for now, inertia in fluid dynamics is considered with a constant
  // density (see if needed / to be debugged)
  multiphysics->set_solution_m1(PhysicsID::VOF, &previous_solutions[0]);

  mass_matrix.reinit(locally_owned_dofs,
                     locally_owned_dofs,
                     dsp,
                     mpi_communicator);

  assemble_mass_matrix_diagonal(mass_matrix);
}

template <int dim>
void
VolumeOfFluid<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*this->fs_mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->VOF,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;

  finish_time_step();
}

template <int dim>
void
VolumeOfFluid<dim>::solve_linear_system(const bool initial_step,
                                        const bool /*renewed_matrix*/)
{
  auto mpi_communicator = triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const double ilu_fill = simulation_parameters.linear_solver.ilu_precond_fill;
  const double ilu_atol = simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol = simulation_parameters.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs, mpi_communicator);

  SolverControl solver_control(
    simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false, simulation_parameters.linear_solver.max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::update_solution_and_constraints(
  TrilinosWrappers::MPI::Vector &solution)
{
  const double penalty_parameter = 100;

  TrilinosWrappers::MPI::Vector lambda(locally_owned_dofs);

  nodal_phase_fraction_owned = solution;

  complete_system_matrix_phase_fraction.residual(lambda,
                                                 nodal_phase_fraction_owned,
                                                 system_rhs_phase_fraction);

  bounding_constraints.clear();
  active_set.clear();

  std::vector<bool> dof_touched(dof_handler.n_dofs(), false);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              Assert(dof_handler.get_fe().dofs_per_cell ==
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcNotImplemented());
              const unsigned int dof_index = cell->vertex_dof_index(v, 0);
              if (locally_owned_dofs.is_element(dof_index))
                {
                  const double solution_value =
                    nodal_phase_fraction_owned(dof_index);
                  if (lambda(dof_index) + penalty_parameter *
                                            mass_matrix(dof_index, dof_index) *
                                            (solution_value - vof_upper_bound) >
                      0)
                    {
                      active_set.add_index(dof_index);
                      bounding_constraints.add_line(dof_index);
                      bounding_constraints.set_inhomogeneity(dof_index,
                                                             vof_upper_bound);
                      nodal_phase_fraction_owned(dof_index) = vof_upper_bound;
                      lambda(dof_index)                     = 0;
                    }
                  else if (lambda(dof_index) +
                             penalty_parameter *
                               mass_matrix(dof_index, dof_index) *
                               (solution_value - vof_lower_bound) <
                           0)
                    {
                      active_set.add_index(dof_index);
                      bounding_constraints.add_line(dof_index);
                      bounding_constraints.set_inhomogeneity(dof_index,
                                                             vof_lower_bound);
                      nodal_phase_fraction_owned(dof_index) = vof_lower_bound;
                      lambda(dof_index)                     = 0;
                    }
                }
            }
        }
    }
  active_set.compress();
  solution = nodal_phase_fraction_owned;
  bounding_constraints.close();
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_L2_projection_interface_sharpening(
  TrilinosWrappers::MPI::Vector &solution)
{
  const double sharpening_threshold =
    this->simulation_parameters.interface_sharpening.sharpening_threshold;
  const double interface_sharpness =
    this->simulation_parameters.interface_sharpening.interface_sharpness;

  FEValues<dim> fe_values_phase_fraction(*this->fs_mapping,
                                         *this->fe,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values |
                                           update_gradients);



  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_q_points    = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_phase_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_phase_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_phase(dofs_per_cell);

  std::vector<double> phase_values(n_q_points);

  system_rhs_phase_fraction    = 0;
  system_matrix_phase_fraction = 0;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_phase_fraction.reinit(cell);

          local_matrix_phase_fraction = 0;
          local_rhs_phase_fraction    = 0;

          fe_values_phase_fraction.get_function_values(solution, phase_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              auto phase_value = phase_values[q];

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_phase[k] = fe_values_phase_fraction.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_phase_fraction(i, j) +=
                        (phi_phase[j] * phi_phase[i]) *
                        fe_values_phase_fraction.JxW(q);
                    }

                  // $$ (if 0 <= \phi <= c)  {\Phi = c ^ (1 - \alpha) * (\phi ^ \alpha)}$$
                  // $$ (if c <  \phi <= 1)  {\Phi = 1 - (1 - c) ^ (1 - \alpha) * (1 - \phi) ^ \alpha}
                  if (phase_value >= 0.0 && phase_value <= sharpening_threshold)
                    local_rhs_phase_fraction(i) +=
                      std::pow(sharpening_threshold,
                               (1 - interface_sharpness)) *
                      std::pow(phase_value, interface_sharpness) *
                      phi_phase[i] * fe_values_phase_fraction.JxW(q);
                  else if (phase_value > sharpening_threshold &&
                           phase_value <= 1.0)
                    {
                      local_rhs_phase_fraction(i) +=
                        (1 -
                         std::pow((1 - sharpening_threshold),
                                  (1 - interface_sharpness)) *
                           std::pow((1 - phase_value), interface_sharpness)) *
                        phi_phase[i] * fe_values_phase_fraction.JxW(q);
                    }
                }
            }
          cell->get_dof_indices(local_dof_indices);
          nonzero_constraints.distribute_local_to_global(
            local_matrix_phase_fraction,
            local_rhs_phase_fraction,
            local_dof_indices,
            system_matrix_phase_fraction,
            system_rhs_phase_fraction);
        }
    }

  system_matrix_phase_fraction.compress(VectorOperation::add);
  system_rhs_phase_fraction.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::solve_interface_sharpening(
  TrilinosWrappers::MPI::Vector &solution)
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-15;

  if (this->simulation_parameters.interface_sharpening.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  TrilinosWrappers::MPI::Vector completely_distributed_phase_fraction_solution(
    this->locally_owned_dofs, triangulation->get_communicator());


  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill =
    this->simulation_parameters.linear_solver.ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_phase_fraction,
                                 preconditionerOptions);

  solver.solve(system_matrix_phase_fraction,
               completely_distributed_phase_fraction_solution,
               system_rhs_phase_fraction,
               *ilu_preconditioner);

  if (this->simulation_parameters.interface_sharpening.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_mass_matrix_diagonal(
  TrilinosWrappers::SparseMatrix &mass_matrix)
{
  QGauss<dim> quadrature_formula(this->cell_quadrature->size());

  FEValues<dim> fe_values(*fe,
                          quadrature_formula,
                          update_values | update_JxW_values);


  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_qpoints     = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_qpoints; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_matrix(i, i) +=
                (fe_values.shape_value(i, q) * fe_values.shape_value(i, q) *
                 fe_values.JxW(q));
          cell->get_dof_indices(local_dof_indices);
          nonzero_constraints.distribute_local_to_global(cell_matrix,
                                                         local_dof_indices,
                                                         mass_matrix);
        }
    }
}

template class VolumeOfFluid<2>;
template class VolumeOfFluid<3>;
