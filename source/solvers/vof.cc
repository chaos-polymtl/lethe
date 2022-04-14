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

#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>
#include <solvers/vof.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>

#include <cmath>

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

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->mapping,
                        dof_handler_fd->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &VolumeOfFluid::assemble_local_system_matrix,
                  &VolumeOfFluid::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_matrix.compress(VectorOperation::add);
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

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_fd);

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
                                              this->system_matrix);
}


template <int dim>
void
VolumeOfFluid<dim>::assemble_system_rhs()
{
  // TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data =
    VOFScratchData<dim>(this->simulation_parameters.physical_properties_manager,
                        *this->fe,
                        *this->cell_quadrature,
                        *this->mapping,
                        dof_handler_fd->get_fe());

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

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_fd);

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
                                              this->system_rhs);
}



template <int dim>
void
VolumeOfFluid<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  data_out.add_data_vector(this->dof_handler, this->present_solution, "phase");
  if (this->simulation_parameters.multiphysics.peeling_wetting)
    {
      // Peeling/wetting output
      data_out.add_data_vector(this->dof_handler, this->marker_pw, "marker_pw");
    }
  if (this->simulation_parameters.multiphysics.continuum_surface_force &&
      simulation_parameters.surface_tension_force.output_VOF_auxiliary_fields)
    {
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        pfg_component_interpretation(
          dim, DataComponentInterpretation::component_is_scalar);
      for (unsigned int i = 0; i < dim; ++i)
        pfg_component_interpretation[i] =
          DataComponentInterpretation::component_is_part_of_vector;

      std::vector<std::string> solution_names(dim, "phase fraction gradient");

      data_out.add_data_vector(pfg_dof_handler,
                               present_pfg_solution,
                               solution_names,
                               pfg_component_interpretation);

      data_out.add_data_vector(curvature_dof_handler,
                               present_curvature_solution,
                               "curvature");
    }
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_L2_error()
{
  auto mpi_communicator = this->triangulation->get_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->error_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  //  Local connectivity
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  const unsigned int  n_q_points = this->error_quadrature->size();
  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->phase;
  exact_solution.set_time(this->simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          fe_values_vof.get_function_values(this->present_solution,
                                            q_scalar_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values_vof.get_quadrature_points(),
                                    q_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = q_scalar_values[q];
              double exact = q_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values_vof.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_volume(int id_fluid_monitored)
{
  auto mpi_communicator = this->triangulation->get_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->error_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  const unsigned int  n_q_points = this->error_quadrature->size();
  std::vector<double> q_scalar_values(n_q_points);

  double volume = 0;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          fe_values_vof.get_function_values(this->present_solution,
                                            q_scalar_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (id_fluid_monitored)
                {
                  case 0:
                    {
                      if (q_scalar_values[q] < 0.5)
                        volume +=
                          fe_values_vof.JxW(q) * (1 - q_scalar_values[q]);
                      break;
                    }
                  case 1:
                    {
                      if (q_scalar_values[q] > 0.5)
                        volume += fe_values_vof.JxW(q) * q_scalar_values[q];
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
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
  auto         mpi_communicator = this->triangulation->get_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        this->error_table.omit_column_from_convergence_rate_evaluation("cells");
      else
        this->error_table.omit_column_from_convergence_rate_evaluation("time");

      this->error_table.set_scientific("error_phase", true);
      this->error_table.set_precision(
        "error_phase", this->simulation_control->get_log_precision());
      this->error_table.write_text(std::cout);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::percolate_time_vectors()
{
  for (unsigned int i = this->previous_solutions.size() - 1; i > 0; --i)
    {
      this->previous_solutions[i] = this->previous_solutions[i - 1];
    }
  this->previous_solutions[0] = this->present_solution;
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

      if (this->simulation_control->is_steady())
        {
          this->error_table.add_value(
            "cells", this->triangulation->n_global_active_cells());
        }
      else
        {
          this->error_table.add_value(
            "time", this->simulation_control->get_current_time());
        }
      this->error_table.add_value("error_phase", phase_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error phase : " << phase_error << std::endl;
        }
    }

  if (simulation_parameters.multiphysics.conservation_monitoring)
    {
      double volume =
        calculate_volume(simulation_parameters.multiphysics.id_fluid_monitored);

      auto         mpi_communicator = this->triangulation->get_communicator();
      unsigned int this_mpi_process(
        Utilities::MPI::this_mpi_process(mpi_communicator));

      if (this_mpi_process == 0)
        {
          // Set conservation monitoring table
          if (this->simulation_control->is_steady())
            {
              this->table_monitoring_vof.add_value(
                "cells", this->triangulation->n_global_active_cells());
            }
          else
            {
              this->table_monitoring_vof.add_value(
                "time", this->simulation_control->get_current_time());
            }

          std::string fluid_id =
            "fluid_" +
            Utilities::int_to_string(
              this->simulation_parameters.multiphysics.id_fluid_monitored, 1);

          this->table_monitoring_vof.add_value("volume_" + fluid_id, volume);
          this->table_monitoring_vof.set_scientific("volume_" + fluid_id, true);

          // Save table to .dat
          std::string filename =
            this->simulation_parameters.simulation_control.output_folder +
            "VOF_monitoring_" + fluid_id + ".dat";
          std::ofstream output(filename.c_str());
          this->table_monitoring_vof.write_text(output);
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::modify_solution()
{
  // Peeling/wetting
  if (this->simulation_parameters.multiphysics.peeling_wetting)
    {
      handle_peeling_wetting();
    }
  // Interface sharpening
  if (this->simulation_parameters.multiphysics.interface_sharpening)
    sharpen_interface();

  if (this->simulation_parameters.multiphysics.continuum_surface_force)
    {
      find_filtered_pfg();
      find_filtered_interface_curvature();
    }
}

template <int dim>
void
VolumeOfFluid<dim>::sharpen_interface()
{
  // Limit the phase fractions between 0 and 1
  update_solution_and_constraints(present_solution);
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    update_solution_and_constraints(previous_solutions[p]);

  // Interface sharpening is done at a constant frequency
  if (this->simulation_control->get_step_number() %
        this->simulation_parameters.interface_sharpening.sharpening_frequency ==
      0)
    {
      if (simulation_parameters.non_linear_solver.verbosity ==
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
            assemble_L2_projection_interface_sharpening(previous_solutions[p]);
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

template <int dim>
void
VolumeOfFluid<dim>::find_filtered_pfg()
{
  assemble_pfg_matrix_and_rhs(present_solution);
  solve_pfg();
}

template <int dim>
void
VolumeOfFluid<dim>::find_filtered_interface_curvature()
{
  assemble_curvature_matrix_and_rhs(present_pfg_solution);
  solve_curvature();
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_pfg_matrix_and_rhs(
  TrilinosWrappers::MPI::Vector &solution)
{
  // Get fe values of VOF phase fraction and phase fraction gradient (pfg)
  FEValues<dim> fe_values_phase_fraction(*this->mapping,
                                         *this->fe,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values |
                                           update_gradients);

  FEValues<dim> fe_values_pfg(*this->mapping,
                              *this->fe_pfg,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values | update_gradients);


  // const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int dofs_per_cell = this->fe_pfg->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_pfg(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_pfg(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<Tensor<1, dim>> phi_pfg(dofs_per_cell);
  std::vector<Tensor<2, dim>> phi_pfg_gradient(dofs_per_cell);

  const FEValuesExtractors::Vector pfg(0);

  std::vector<Tensor<1, dim>> phase_gradient_values(n_q_points);
  std::vector<Tensor<1, dim>> pfg_values(n_q_points);
  std::vector<Tensor<2, dim>> pfg_gradient_values(n_q_points);

  // Reinitialize system matrix and rhs for the pfg
  system_rhs_pfg    = 0;
  system_matrix_pfg = 0;

  for (const auto &pfg_cell : this->pfg_dof_handler.active_cell_iterators())
    {
      if (pfg_cell->is_locally_owned())
        {
          // Gather the active cell iterator related to the VOF phase fraction
          typename DoFHandler<dim>::active_cell_iterator cell(
            &(*this->triangulation),
            pfg_cell->level(),
            pfg_cell->index(),
            &this->dof_handler);

          fe_values_phase_fraction.reinit(cell);
          fe_values_pfg.reinit(pfg_cell);

          local_matrix_pfg = 0;
          local_rhs_pfg    = 0;

          // Get phase fraction values, pfg values and gradients
          fe_values_phase_fraction.get_function_gradients(
            solution, phase_gradient_values);

          fe_values_pfg[pfg].get_function_values(present_pfg_solution,
                                                 pfg_values);

          fe_values_pfg[pfg].get_function_gradients(present_pfg_solution,
                                                    pfg_gradient_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_pfg[k]          = fe_values_pfg[pfg].value(k, q);
                  phi_pfg_gradient[k] = fe_values_pfg[pfg].gradient(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_pfg(i, j) +=
                        (phi_pfg[j] * phi_pfg[i] +
                         simulation_parameters.surface_tension_force
                             .phase_fraction_gradient_filter_value *
                           scalar_product(phi_pfg_gradient[i],
                                          phi_pfg_gradient[j])) *
                        fe_values_pfg.JxW(q);
                    }

                  // rhs
                  local_rhs_pfg(i) -=
                    (phi_pfg[i] * pfg_values[q] +
                     simulation_parameters.surface_tension_force
                         .phase_fraction_gradient_filter_value *
                       scalar_product(phi_pfg_gradient[i],
                                      pfg_gradient_values[q]) -
                     phi_pfg[i] * phase_gradient_values[q]) *
                    fe_values_pfg.JxW(q);
                }
            }

          pfg_cell->get_dof_indices(local_dof_indices);
          pfg_constraints.distribute_local_to_global(local_matrix_pfg,
                                                     local_rhs_pfg,
                                                     local_dof_indices,
                                                     system_matrix_pfg,
                                                     system_rhs_pfg);
        }
    }
  system_matrix_pfg.compress(VectorOperation::add);
  system_rhs_pfg.compress(VectorOperation::add);
}


template <int dim>
void
VolumeOfFluid<dim>::solve_pfg()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-10;

  TrilinosWrappers::MPI::Vector completely_distributed_pfg_solution(
    this->locally_owned_dofs_pfg, triangulation->get_communicator());

  completely_distributed_pfg_solution = present_pfg_solution;

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  const double ilu_fill =
    this->simulation_parameters.linear_solver.ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_pfg, preconditionerOptions);
  solver.solve(system_matrix_pfg,
               completely_distributed_pfg_solution,
               system_rhs_pfg,
               *ilu_preconditioner);

  if (this->simulation_parameters.surface_tension_force.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver (phase fraction gradient) took : "
                  << solver_control.last_step() << " steps " << std::endl;
    }

  pfg_constraints.distribute(completely_distributed_pfg_solution);

  present_pfg_solution = completely_distributed_pfg_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_curvature_matrix_and_rhs(
  TrilinosWrappers::MPI::Vector &present_pfg_solution)
{
  // Get fe values of phase fraction gradient (pfg) and curvature
  FEValues<dim> fe_values_curvature(*this->curvature_mapping,
                                    *this->fe_curvature,
                                    *this->cell_quadrature,
                                    update_values | update_quadrature_points |
                                      update_JxW_values | update_gradients);

  FEValues<dim> fe_values_pfg(*this->mapping,
                              *this->fe_pfg,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_curvature->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_curvature(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_curvature(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double>         phi_curvature(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_curvature_gradient(dofs_per_cell);

  pfg_values       = std::vector<Tensor<1, dim>>(n_q_points);
  curvature_values = std::vector<double>(n_q_points);
  std::vector<Tensor<1, dim>>      curvature_gradient_values(n_q_points);
  const FEValuesExtractors::Vector pfg(0);

  // Reinitialize system matrix and rhs for the curvature
  system_rhs_curvature    = 0;
  system_matrix_curvature = 0;

  for (const auto &curvature_cell :
       this->curvature_dof_handler.active_cell_iterators())
    {
      if (curvature_cell->is_locally_owned())
        {
          // Gather the active cell iterator related to the phase fraction
          // gradient (pfg)
          typename DoFHandler<dim>::active_cell_iterator pfg_cell(
            &(*this->triangulation),
            curvature_cell->level(),
            curvature_cell->index(),
            &this->pfg_dof_handler);

          fe_values_pfg.reinit(pfg_cell);
          fe_values_curvature.reinit(curvature_cell);

          local_matrix_curvature = 0;
          local_rhs_curvature    = 0;

          // Get pfg values, curvature values and gradients
          fe_values_pfg[pfg].get_function_values(present_pfg_solution,
                                                 pfg_values);

          fe_values_curvature.get_function_values(present_curvature_solution,
                                                  curvature_values);

          fe_values_curvature.get_function_gradients(present_curvature_solution,
                                                     curvature_gradient_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_curvature[k] = fe_values_curvature.shape_value(k, q);
                  phi_curvature_gradient[k] =
                    fe_values_curvature.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_curvature(i, j) +=
                        (phi_curvature[j] * phi_curvature[i] +
                         simulation_parameters.surface_tension_force
                             .curvature_filter_value *
                           scalar_product(phi_curvature_gradient[i],
                                          phi_curvature_gradient[j])) *
                        fe_values_curvature.JxW(q);
                    }
                  // rhs
                  local_rhs_curvature(i) -=
                    (phi_curvature[i] * curvature_values[q] +
                     simulation_parameters.surface_tension_force
                         .curvature_filter_value *
                       scalar_product(phi_curvature_gradient[i],
                                      curvature_gradient_values[q]) -
                     phi_curvature_gradient[i] *
                       (pfg_values[q] / (pfg_values[q].norm() + DBL_MIN))) *
                    fe_values_curvature.JxW(q);
                }
            }

          curvature_cell->get_dof_indices(local_dof_indices);
          curvature_constraints.distribute_local_to_global(
            local_matrix_curvature,
            local_rhs_curvature,
            local_dof_indices,
            system_matrix_curvature,
            system_rhs_curvature);
        }
    }
  system_matrix_curvature.compress(VectorOperation::add);
  system_rhs_curvature.compress(VectorOperation::add);
}

template <int dim>
void
VolumeOfFluid<dim>::solve_curvature()
{
  const double linear_solver_tolerance = 1e-10;

  TrilinosWrappers::MPI::Vector completely_distributed_curvature_solution(
    this->locally_owned_dofs_curvature, triangulation->get_communicator());

  completely_distributed_curvature_solution = present_curvature_solution;

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverCG solver(solver_control);

  const double ilu_fill =
    this->simulation_parameters.linear_solver.ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.ilu_precond_rtol;

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix_curvature,
                                 preconditionerOptions);
  solver.solve(system_matrix_curvature,
               completely_distributed_curvature_solution,
               system_rhs_curvature,
               *ilu_preconditioner);

  if (this->simulation_parameters.surface_tension_force.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << " -Iterative solver (curvature) took : "
                  << solver_control.last_step() << " steps " << std::endl;
    }

  curvature_constraints.distribute(completely_distributed_curvature_solution);

  present_curvature_solution = completely_distributed_curvature_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::pre_mesh_adaptation()
{
  this->solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);

  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      this->previous_solutions_transfer[i]
        .prepare_for_coarsening_and_refinement(this->previous_solutions[i]);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = this->triangulation->get_communicator();

  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::Vector tmp(this->locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  this->solution_transfer.interpolate(tmp);

  // Distribute constraints
  this->nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  this->present_solution = tmp;

  // Transfer previous solutions
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      TrilinosWrappers::MPI::Vector tmp_previous_solution(
        this->locally_owned_dofs, mpi_communicator);
      this->previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      this->nonzero_constraints.distribute(tmp_previous_solution);
      this->previous_solutions[i] = tmp_previous_solution;
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
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(phase));
    }
}

template <int dim>
void
VolumeOfFluid<dim>::write_checkpoint()
{
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;

  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }
  this->solution_transfer.prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
VolumeOfFluid<dim>::read_checkpoint()
{
  auto mpi_communicator        = this->triangulation->get_communicator();
  auto previous_solutions_size = this->previous_solutions.size();
  this->pcout << "Reading VOF checkpoint" << std::endl;

  std::vector<TrilinosWrappers::MPI::Vector *> input_vectors(
    1 + previous_solutions_size);
  TrilinosWrappers::MPI::Vector distributed_system(this->locally_owned_dofs,
                                                   mpi_communicator);
  input_vectors[0] = &distributed_system;


  std::vector<TrilinosWrappers::MPI::Vector> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions_size);
  for (unsigned int i = 0; i < previous_solutions_size; ++i)
    {
      distributed_previous_solutions.emplace_back(
        TrilinosWrappers::MPI::Vector(this->locally_owned_dofs,
                                      mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  this->solution_transfer.deserialize(input_vectors);

  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions_size; ++i)
    {
      this->previous_solutions[i] = distributed_previous_solutions[i];
    }
}


template <int dim>
void
VolumeOfFluid<dim>::setup_dofs()
{
  auto mpi_communicator = triangulation->get_communicator();

  pfg_dof_handler.distribute_dofs(*fe_pfg);

  locally_owned_dofs_pfg = pfg_dof_handler.locally_owned_dofs();

  DoFTools::extract_locally_relevant_dofs(pfg_dof_handler,
                                          locally_relevant_dofs_pfg);

  pfg_constraints.clear();
  pfg_constraints.reinit(locally_relevant_dofs_pfg);
  DoFTools::make_hanging_node_constraints(pfg_dof_handler, pfg_constraints);
  pfg_constraints.close();

  nodal_pfg_relevant.reinit(locally_owned_dofs_pfg,
                            locally_relevant_dofs_pfg,
                            mpi_communicator);

  nodal_pfg_owned.reinit(locally_owned_dofs_pfg, mpi_communicator);



  DynamicSparsityPattern dsp_pfg(locally_relevant_dofs_pfg);
  DoFTools::make_sparsity_pattern(pfg_dof_handler,
                                  dsp_pfg,
                                  pfg_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(dsp_pfg,
                                             locally_owned_dofs_pfg,
                                             mpi_communicator,
                                             locally_relevant_dofs_pfg);


  // Initialization of phase fraction gradient matrice and rhs for the
  // calculation surface tension force
  system_matrix_pfg.reinit(locally_owned_dofs_pfg,
                           locally_owned_dofs_pfg,
                           dsp_pfg,
                           mpi_communicator);

  system_rhs_pfg.reinit(locally_owned_dofs_pfg, mpi_communicator);

  present_pfg_solution.reinit(locally_owned_dofs_pfg,
                              locally_relevant_dofs_pfg,
                              mpi_communicator);

  // Curvature
  curvature_dof_handler.distribute_dofs(*fe_curvature);

  locally_owned_dofs_curvature = curvature_dof_handler.locally_owned_dofs();

  DoFTools::extract_locally_relevant_dofs(curvature_dof_handler,
                                          locally_relevant_dofs_curvature);

  curvature_constraints.clear();
  curvature_constraints.reinit(locally_relevant_dofs_curvature);
  DoFTools::make_hanging_node_constraints(curvature_dof_handler,
                                          curvature_constraints);
  curvature_constraints.close();

  nodal_curvature_relevant.reinit(locally_owned_dofs_curvature,
                                  locally_relevant_dofs_curvature,
                                  mpi_communicator);

  nodal_curvature_owned.reinit(locally_owned_dofs_curvature, mpi_communicator);



  DynamicSparsityPattern dsp_curvature(locally_relevant_dofs_curvature);
  DoFTools::make_sparsity_pattern(curvature_dof_handler,
                                  dsp_curvature,
                                  curvature_constraints,
                                  false);

  SparsityTools::distribute_sparsity_pattern(dsp_curvature,
                                             locally_owned_dofs_curvature,
                                             mpi_communicator,
                                             locally_relevant_dofs_curvature);


  // Initialization of curvature matrice and rhs for the
  // calculation surface tension force
  system_matrix_curvature.reinit(locally_owned_dofs_curvature,
                                 locally_owned_dofs_curvature,
                                 dsp_curvature,
                                 mpi_communicator);

  system_rhs_curvature.reinit(locally_owned_dofs_curvature, mpi_communicator);

  present_curvature_solution.reinit(locally_owned_dofs_curvature,
                                    locally_relevant_dofs_curvature,
                                    mpi_communicator);

  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                mpi_communicator);

  // Previous solutions for transient schemes
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      mpi_communicator);
    }

  this->system_rhs.reinit(this->locally_owned_dofs, mpi_communicator);

  this->newton_update.reinit(this->locally_owned_dofs, mpi_communicator);

  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      mpi_communicator);

  {
    this->nonzero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->nonzero_constraints);
  }
  this->nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    this->zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            this->zero_constraints);
  }
  this->zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(this->dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  this->nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             this->locally_owned_dofs,
                                             mpi_communicator,
                                             this->locally_relevant_dofs);

  // Initialization of phase fraction matrices for interface sharpening.
  // system_matrix_phase_fraction is used in
  // assemble_L2_projection_interface_sharpening for assembling the system for
  // sharpening the interface, while complete_system_matrix_phase_fraction is
  // used in update_solution_and_constraints to limit the phase fraction values
  // between 0 and 1. Accoring to step-41, to limit the phase fractions we
  // compute the Lagrange multiplier as the residual of the original linear
  // system, given via the variables complete_system_matrix_phase_fraction and
  // complete_system_rhs_phase_fraction
  system_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                      this->locally_owned_dofs,
                                      dsp,
                                      mpi_communicator);

  complete_system_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                               this->locally_owned_dofs,
                                               dsp,
                                               mpi_communicator);

  complete_system_rhs_phase_fraction.reinit(this->locally_owned_dofs,
                                            mpi_communicator);

  // In update_solution_and_constraints (which limits the phase fraction
  // between 0 and 1) nodal_phase_fraction_owned copies the solution, then
  // limits it, and finally updates (rewrites) the solution.
  nodal_phase_fraction_owned.reinit(this->locally_owned_dofs, mpi_communicator);

  // Right hand side of the interface sharpening problem (used in
  // assemble_L2_projection_interface_sharpening).
  system_rhs_phase_fraction.reinit(this->locally_owned_dofs, mpi_communicator);

  this->system_matrix.reinit(this->locally_owned_dofs,
                             this->locally_owned_dofs,
                             dsp,
                             mpi_communicator);

  // Initialize peeling/wetting variables
  marker_pw.reinit(locally_owned_dofs, mpi_communicator);
  dofs_wet.reinit(locally_owned_dofs, mpi_communicator);
  dofs_peeled.reinit(locally_owned_dofs, mpi_communicator);

  this->pcout << "   Number of VOF degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;

  // Provide the VOF dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::VOF, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::VOF, &this->present_solution);

  // the fluid at present iteration is solved BEFORE the VOF (see map
  // solve_pre_fluid defined in multiphysics_interface.h), and after percolate
  // is called for the previous iteration.
  // NB: for now, inertia in fluid dynamics is considered with a constant
  // density (see if needed / to be debugged)
  multiphysics->set_solution_m1(PhysicsID::VOF, &this->previous_solutions[0]);

  mass_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                    this->locally_owned_dofs,
                                    dsp,
                                    mpi_communicator);

  assemble_mass_matrix_diagonal(mass_matrix_phase_fraction);
}

template <int dim>
void
VolumeOfFluid<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*this->mapping,
                           this->dof_handler,
                           simulation_parameters.initial_condition->VOF,
                           this->newton_update);
  this->nonzero_constraints.distribute(this->newton_update);
  this->present_solution = this->newton_update;

  finish_time_step();
}

template <int dim>
void
VolumeOfFluid<dim>::solve_linear_system(const bool initial_step,
                                        const bool /*renewed_matrix*/)
{
  auto mpi_communicator = this->triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? this->nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * this->system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.non_linear_solver.verbosity !=
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

  ilu_preconditioner.initialize(this->system_matrix, preconditionerOptions);

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, mpi_communicator);

  SolverControl solver_control(
    simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false, simulation_parameters.linear_solver.max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.non_linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

// This function is explained in detail in step-41 of deal.II tutorials
template <int dim>
void
VolumeOfFluid<dim>::update_solution_and_constraints(
  TrilinosWrappers::MPI::Vector &solution)
{
  // This is a penalty parameter for limiting the phase fraction
  // in the range of [0,1]. According to step 41, this parameter depends
  // on the problem itself and needs to be chosen large enough (for example
  // there is no convergence using the penalty_parameter = 1)
  const double penalty_parameter = 100;

  TrilinosWrappers::MPI::Vector lambda(this->locally_owned_dofs);

  nodal_phase_fraction_owned = solution;

  complete_system_matrix_phase_fraction.residual(lambda,
                                                 nodal_phase_fraction_owned,
                                                 system_rhs_phase_fraction);

  this->bounding_constraints.clear();

  std::vector<bool> dof_touched(this->dof_handler.n_dofs(), false);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            {
              Assert(this->dof_handler.get_fe().dofs_per_cell ==
                       GeometryInfo<dim>::vertices_per_cell,
                     ExcNotImplemented());
              const unsigned int dof_index = cell->vertex_dof_index(v, 0);
              if (this->locally_owned_dofs.is_element(dof_index))
                {
                  const double solution_value =
                    nodal_phase_fraction_owned(dof_index);
                  if (lambda(dof_index) +
                        penalty_parameter *
                          mass_matrix_phase_fraction(dof_index, dof_index) *
                          (solution_value - this->phase_upper_bound) >
                      0)
                    {
                      this->bounding_constraints.add_line(dof_index);
                      this->bounding_constraints.set_inhomogeneity(
                        dof_index, this->phase_upper_bound);
                      nodal_phase_fraction_owned(dof_index) =
                        this->phase_upper_bound;
                      lambda(dof_index) = 0;
                    }
                  else if (lambda(dof_index) +
                             penalty_parameter *
                               mass_matrix_phase_fraction(dof_index,
                                                          dof_index) *
                               (solution_value - this->phase_lower_bound) <
                           0)
                    {
                      this->bounding_constraints.add_line(dof_index);
                      this->bounding_constraints.set_inhomogeneity(
                        dof_index, this->phase_lower_bound);
                      nodal_phase_fraction_owned(dof_index) =
                        this->phase_lower_bound;
                      lambda(dof_index) = 0;
                    }
                }
            }
        }
    }
  solution = nodal_phase_fraction_owned;
  this->bounding_constraints.close();
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

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values | update_gradients);



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
          fe_values_vof.reinit(cell);

          local_matrix_phase_fraction = 0;
          local_rhs_phase_fraction    = 0;

          fe_values_vof.get_function_values(solution, phase_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              auto phase_value = phase_values[q];

              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_phase[k] = fe_values_vof.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_phase_fraction(i, j) +=
                        (phi_phase[j] * phi_phase[i]) * fe_values_vof.JxW(q);
                    }

                  // $$ (if 0 <= \phi <= c)  {\Phi = c ^ (1 - \alpha) * (\phi ^
                  // \alpha)}$$
                  // $$ (if c <  \phi <= 1)  {\Phi = 1 - (1 - c) ^ (1 - \alpha)
                  // * (1 - \phi) ^ \alpha}
                  if (phase_value >= 0.0 && phase_value <= sharpening_threshold)
                    local_rhs_phase_fraction(i) +=
                      std::pow(sharpening_threshold,
                               (1 - interface_sharpness)) *
                      std::pow(phase_value, interface_sharpness) *
                      phi_phase[i] * fe_values_vof.JxW(q);
                  else
                    {
                      local_rhs_phase_fraction(i) +=
                        (1 -
                         std::pow((1 - sharpening_threshold),
                                  (1 - interface_sharpness)) *
                           std::pow((1 - phase_value), interface_sharpness)) *
                        phi_phase[i] * fe_values_vof.JxW(q);
                    }
                }
            }
          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
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

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

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

  this->nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}

// This function is explained in detail in step-41 of deal.II tutorials:
// We get the mass matrix to be diagonal by choosing the trapezoidal rule
// for quadrature. Doing so we do not really need the triple loop over
// quadrature points, indices i and indices j any more and can, instead, just
// use a double loop.
template <int dim>
void
VolumeOfFluid<dim>::assemble_mass_matrix_diagonal(
  TrilinosWrappers::SparseMatrix &mass_matrix)
{
  QGauss<dim> quadrature_formula(this->cell_quadrature->size());

  FEValues<dim> fe_values_vof(*this->fe,
                              quadrature_formula,
                              update_values | update_JxW_values);


  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_qpoints     = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_qpoints; ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              cell_matrix(i, i) +=
                (fe_values_vof.shape_value(i, q) *
                 fe_values_vof.shape_value(i, q) * fe_values_vof.JxW(q));
          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            cell_matrix, local_dof_indices, mass_matrix);
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::handle_peeling_wetting()
{
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_vof.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions_vof.type[i_bc] ==
          BoundaryConditions::BoundaryType::pw)
        {
          // Parse fluid present solution to apply_peeling_wetting method
          if (multiphysics->fluid_dynamics_is_block())
            {
              const TrilinosWrappers::MPI::BlockVector current_solution_fd(
                *multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
              // apply_peeling_wetting is templated with
              // current_solution_fd VectorType
              apply_peeling_wetting(i_bc, current_solution_fd);
            }
          else
            {
              const TrilinosWrappers::MPI::Vector current_solution_fd(
                *multiphysics->get_solution(PhysicsID::fluid_dynamics));
              // apply_peeling_wetting is templated with
              // current_solution_fd VectorType
              apply_peeling_wetting(i_bc, current_solution_fd);
            }
        }
    } // end loop on boundary_conditions_vof

  // Output total of peeled/wet cells in the entire domain
  if (this->simulation_parameters.non_linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      int nb_cells_wet = this->dofs_wet.l1_norm() / this->fe->dofs_per_cell;
      int nb_cells_peeled =
        this->dofs_peeled.l1_norm() / this->fe->dofs_per_cell;

      this->pcout << "Peeling/wetting correction at step "
                  << this->simulation_control->get_step_number() << std::endl;
      //      this->pcout << "  -number of wet cells: " << this->nb_cells_wet
      this->pcout << "  -number of wet cells: " << nb_cells_wet << std::endl;
      //      this->pcout << "  -number of peeled cells: " <<
      //      this->nb_cells_peeled
      this->pcout << "  -number of peeled cells: " << nb_cells_peeled
                  << std::endl;
    }
}

template <int dim>
template <typename VectorType>
void
VolumeOfFluid<dim>::apply_peeling_wetting(const unsigned int i_bc,
                                          const VectorType &current_solution_fd)
{
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  // Initializations
  locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  std::vector<types::global_dof_index> dof_indices_vof(
    fe->dofs_per_cell); //  local connectivity
  auto mpi_communicator = this->triangulation->get_communicator();

  solution_pw.reinit(this->locally_owned_dofs, mpi_communicator);
  solution_pw = this->present_solution;

  FEFaceValues<dim> fe_face_values_vof(*this->mapping,
                                       *this->fe,
                                       *this->face_quadrature,
                                       update_values |
                                         update_quadrature_points);

  FEFaceValues<dim> fe_face_values_fd(*this->mapping,
                                      dof_handler_fd->get_fe(),
                                      *this->face_quadrature,
                                      update_values | update_quadrature_points |
                                        update_gradients);

  const unsigned int n_q_points = this->cell_quadrature->size();

  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>      pressure_gradients(n_q_points);
  std::vector<double>              phase_values(n_q_points);

  unsigned int boundary_id =
    this->simulation_parameters.boundary_conditions.id[i_bc];

  // Physical properties
  const auto density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::map<field, std::vector<double>> fields;

  std::vector<double> density_0(n_q_points);
  std::vector<double> density_1(n_q_points);

  // Useful definitions for readability
  const double wetting_threshold =
    this->simulation_parameters.boundary_conditions_vof.wetting_threshold[i_bc];
  const double peeling_threshold =
    this->simulation_parameters.boundary_conditions_vof.peeling_threshold[i_bc];

  // Loop on cell_vof
  for (const auto &cell_vof : dof_handler.active_cell_iterators())
    {
      if (cell_vof->is_locally_owned() && cell_vof->at_boundary())
        {
          // Initialize average values on the cell
          double       id_denser_fluid_cell(0);
          double       id_lighter_fluid_cell(0);
          double       phase_values_cell(0);
          double       pressure_values_cell(0);
          unsigned int nb_pressure_grad_meet_peel_condition(0);

          // Local index (see deal.II step 4)
          cell_vof->get_dof_indices(dof_indices_vof);

          for (const auto face : cell_vof->face_indices())
            { // TODO ajouter condition que soit sur la BC de pw
              if (cell_vof->face(face)->at_boundary() &&
                  cell_vof->face(face)->boundary_id() == boundary_id)
                {
                  // Get fluid dynamics active cell iterator
                  typename DoFHandler<dim>::active_cell_iterator cell_fd(
                    &(*(this->triangulation)),
                    cell_vof->level(),
                    cell_vof->index(),
                    dof_handler_fd);

                  // Reinit fe_face_values for Fluid Dynamics and VOF
                  fe_face_values_vof.reinit(cell_vof, face);
                  fe_face_values_fd.reinit(cell_fd, face);

                  // Get pressure (values, gradient)
                  fe_face_values_fd[pressure].get_function_values(
                    current_solution_fd, pressure_values);
                  fe_face_values_fd[pressure].get_function_gradients(
                    current_solution_fd, pressure_gradients);

                  // Get phase values
                  fe_face_values_vof.get_function_values(present_solution,
                                                         phase_values);

                  // Calculate physical properties for the cell
                  density_models[0]->vector_value(fields, density_0);
                  density_models[1]->vector_value(fields, density_1);

                  // Loop on the quadrature points
                  for (unsigned int q = 0; q < n_q_points; q++)
                    {
                      // Get denser/lighter fluid id
                      if (density_1[q] > density_0[q])
                        {
                          id_denser_fluid_cell += 1;
                        }
                      else
                        {
                          id_lighter_fluid_cell += 1;
                        }

                      // Peeling condition on the pressure gradient
                      for (unsigned int d = 0; d < dim; d++)
                        {
                          if (pressure_gradients[q][d] < peeling_threshold)
                            {
                              nb_pressure_grad_meet_peel_condition += 1;
                            }
                        }

                      pressure_values_cell += pressure_values[q];
                      phase_values_cell += phase_values[q];

                    } // end loop on quadrature points
                }
            } // end loop on faces

          // Caculate average values on the cell
          id_denser_fluid_cell  = round(id_denser_fluid_cell / n_q_points);
          id_lighter_fluid_cell = round(id_lighter_fluid_cell / n_q_points);
          pressure_values_cell  = pressure_values_cell / n_q_points;
          phase_values_cell     = phase_values_cell / n_q_points;

          // Wetting of lower density fluid
          if ((pressure_values_cell > wetting_threshold) &&
              ((id_denser_fluid_cell == 1 && phase_values_cell < 0.5) ||
               (id_denser_fluid_cell == 0 && phase_values_cell > 0.5)))
            {
              change_cell_phase(PhaseChange::wetting,
                                id_denser_fluid_cell,
                                solution_pw,
                                dof_indices_vof);
            }

          // Peeling of higher density fluid
          else if ((nb_pressure_grad_meet_peel_condition >
                    dim * n_q_points / 2) &&
                   ((id_denser_fluid_cell == 1 && phase_values_cell > 0.5) ||
                    (id_denser_fluid_cell == 0 && phase_values_cell < 0.5)))
            {
              change_cell_phase(PhaseChange::peeling,
                                id_lighter_fluid_cell,
                                solution_pw,
                                dof_indices_vof);
            }

        } // end condition cell at boundary
    }     // end loop on cells

  present_solution = solution_pw;
}

template void
VolumeOfFluid<2>::apply_peeling_wetting<TrilinosWrappers::MPI::Vector>(
  const unsigned int                   i_bc,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
VolumeOfFluid<3>::apply_peeling_wetting<TrilinosWrappers::MPI::Vector>(
  const unsigned int                   i_bc,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
VolumeOfFluid<2>::apply_peeling_wetting<TrilinosWrappers::MPI::BlockVector>(
  const unsigned int                        i_bc,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template void
VolumeOfFluid<3>::apply_peeling_wetting<TrilinosWrappers::MPI::BlockVector>(
  const unsigned int                        i_bc,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template <int dim>
void
VolumeOfFluid<dim>::change_cell_phase(
  const PhaseChange &                         type,
  const unsigned int &                        new_phase,
  TrilinosWrappers::MPI::Vector &             solution_pw,
  const std::vector<types::global_dof_index> &dof_indices_vof)
{
  // TMP
  //  std::cout << "nb_cells_peeled entree change_cell_phase"
  //            << this->nb_cells_peeled << std::endl;

  if (type == PhaseChange::wetting)
    {
      for (unsigned int k = 0; k < fe->dofs_per_cell; ++k)
        {
          this->marker_pw[dof_indices_vof[k]] = 1;
          solution_pw[dof_indices_vof[k]]     = new_phase;
          this->dofs_wet                      = 1;
        }
      //      // increment cells count
      //      this->nb_cells_wet++;
    }
  else if (type == PhaseChange::peeling)
    {
      for (unsigned int k = 0; k < fe->dofs_per_cell; ++k)
        {
          this->marker_pw[dof_indices_vof[k]] = -1;
          solution_pw[dof_indices_vof[k]]     = new_phase;
          this->dofs_peeled                   = 1;
        }
      //      // increment cells count
      //      this->nb_cells_peeled++;

      //      // TMP
      //      std::cout << "nb_cells_peeled post peeling" <<
      //      this->nb_cells_peeled
      //                << std::endl;
    }
}


template class VolumeOfFluid<2>;
template class VolumeOfFluid<3>;
