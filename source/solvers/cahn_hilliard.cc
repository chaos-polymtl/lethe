// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/cahn_hilliard.h>
#include <solvers/cahn_hilliard_assemblers.h>
#include <solvers/cahn_hilliard_scratch_data.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
void
CahnHilliard<dim>::setup_assemblers()
{
  AssertThrow(
    is_sdirk(this->simulation_control->get_assembly_method()) == false,
    ExcMessage("The SDIRK scheme is not yet supported for this physics"));

  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.emplace_back(
        std::make_shared<CahnHilliardAssemblerBDF<dim>>(
          this->simulation_control));
    }

  // For all the assemblers below, the parameter epsilon is passed to the
  // constructor explicitly because it might be calculated directly using
  // properties from the triangulation (e.g. minimal cell-size). Consequently,
  // it cannot be used directly from the simulation parameters and it is passed
  // to the constructor separately.

  // Angle of contact boundary condition
  this->assemblers.emplace_back(
    std::make_shared<CahnHilliardAssemblerAngleOfContact<dim>>(
      this->simulation_control,
      this->simulation_parameters.multiphysics.cahn_hilliard_parameters,
      (this->simulation_parameters.multiphysics.cahn_hilliard_parameters
         .epsilon_set_method == Parameters::EpsilonSetMethod::manual) ?
        this->simulation_parameters.multiphysics.cahn_hilliard_parameters
          .epsilon :
        GridTools::minimal_cell_diameter(*triangulation),
      this->simulation_parameters.boundary_conditions_cahn_hilliard));

  // Free angle of contact boundary condition
  this->assemblers.emplace_back(
    std::make_shared<CahnHilliardAssemblerFreeAngle<dim>>(
      this->simulation_control,
      this->simulation_parameters.multiphysics.cahn_hilliard_parameters,
      (this->simulation_parameters.multiphysics.cahn_hilliard_parameters
         .epsilon_set_method == Parameters::EpsilonSetMethod::manual) ?
        this->simulation_parameters.multiphysics.cahn_hilliard_parameters
          .epsilon :
        GridTools::minimal_cell_diameter(*triangulation),
      this->simulation_parameters.boundary_conditions_cahn_hilliard));

  // Core assembler
  this->assemblers.emplace_back(
    std::make_shared<CahnHilliardAssemblerCore<dim>>(
      this->simulation_control,
      this->simulation_parameters.multiphysics.cahn_hilliard_parameters,
      (this->simulation_parameters.multiphysics.cahn_hilliard_parameters
         .epsilon_set_method == Parameters::EpsilonSetMethod::manual) ?
        this->simulation_parameters.multiphysics.cahn_hilliard_parameters
          .epsilon :
        GridTools::minimal_cell_diameter(*triangulation)));
}

template <int dim>
void
CahnHilliard<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = CahnHilliardScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &CahnHilliard::assemble_local_system_matrix,
                  &CahnHilliard::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
CahnHilliard<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  CahnHilliardScratchData<dim>                         &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!copy_data.cell_is_local)
    return;

  auto source_term = simulation_parameters.source_term.cahn_hilliard_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
    }
  else
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
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
CahnHilliard<dim>::copy_local_matrix_to_global_matrix(
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
CahnHilliard<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = CahnHilliardScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &CahnHilliard::assemble_local_system_rhs,
                  &CahnHilliard::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
CahnHilliard<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  CahnHilliardScratchData<dim>                         &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!copy_data.cell_is_local)
    return;

  auto source_term = simulation_parameters.source_term.cahn_hilliard_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
    }
  else
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
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
CahnHilliard<dim>::copy_local_rhs_to_global_rhs(
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
CahnHilliard<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  // Add the interpretation of the solution. The first component is the
  // phase order (Phi) and the following one is the chemical potential (eta)

  std::vector<std::string> solution_names;
  solution_names.emplace_back("phase_order");
  solution_names.emplace_back("chemical_potential");

  std::vector<std::string> solution_names_filtered;
  solution_names_filtered.emplace_back("phase_order_filtered");
  solution_names_filtered.emplace_back("chemical_potential_filtered");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      2, DataComponentInterpretation::component_is_scalar);

  data_out.add_data_vector(dof_handler,
                           present_solution,
                           solution_names,
                           data_component_interpretation);

  // Filter phase fraction
  data_out.add_data_vector(dof_handler,
                           filtered_solution,
                           solution_names_filtered,
                           data_component_interpretation);
}

template <int dim>
std::pair<double, double>
CahnHilliard<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*this->mapping,
                          *fe,
                          *this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar phase_order(0);
  const FEValuesExtractors::Scalar chemical_potential(1);
  const unsigned int               n_q_points = this->cell_quadrature->size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(2));

  std::vector<double> local_phase_order_values(n_q_points);
  std::vector<double> local_chemical_potential_values(n_q_points);

  auto &exact_solution =
    simulation_parameters.analytical_solution->cahn_hilliard;
  exact_solution.set_time(simulation_control->get_current_time());

  double l2_error_phase_order        = 0.;
  double l2_error_chemical_potential = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[phase_order].get_function_values(evaluation_point,
                                                     local_phase_order_values);
          fe_values[chemical_potential].get_function_values(
            evaluation_point, local_chemical_potential_values);

          // Get the exact solution at all Gauss points
          exact_solution.vector_value_list(fe_values.get_quadrature_points(),
                                           q_exactSol);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Find the values of x and u_h (the finite element solution) at
              // the quadrature points
              double phase_order_sim   = local_phase_order_values[q];
              double phase_order_exact = q_exactSol[q][0];

              l2_error_phase_order += (phase_order_sim - phase_order_exact) *
                                      (phase_order_sim - phase_order_exact) *
                                      fe_values.JxW(q);

              double chemical_potential_sim =
                local_chemical_potential_values[q];
              double chemical_potential_exact = q_exactSol[q][1];

              l2_error_chemical_potential +=
                (chemical_potential_sim - chemical_potential_exact) *
                (chemical_potential_sim - chemical_potential_exact) *
                fe_values.JxW(q);
            }
        }
    }
  l2_error_phase_order =
    Utilities::MPI::sum(l2_error_phase_order, mpi_communicator);
  l2_error_chemical_potential =
    Utilities::MPI::sum(l2_error_chemical_potential, mpi_communicator);

  return std::make_pair(std::sqrt(l2_error_phase_order),
                        std::sqrt(l2_error_chemical_potential));
}


template <int dim>
void
CahnHilliard<dim>::calculate_phase_statistics()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar phase_order(0);

  const unsigned int          n_q_points = cell_quadrature->size();
  std::vector<double>         local_phase_order_values(n_q_points);
  std::vector<Tensor<1, dim>> local_phase_order_gradients(n_q_points);

  double integral(0.);
  double max_phase_value(std::numeric_limits<double>::min());
  double min_phase_value(std::numeric_limits<double>::max());
  double volume_0(0.);
  double volume_1(0.);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[phase_order].get_function_values(present_solution,
                                                     local_phase_order_values);
          fe_values[phase_order].get_function_gradients(
            present_solution, local_phase_order_gradients);
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              integral += local_phase_order_values[q] * fe_values.JxW(q);
              max_phase_value =
                std::max(local_phase_order_values[q], max_phase_value);
              min_phase_value =
                std::min(local_phase_order_values[q], min_phase_value);
              volume_0 +=
                (1 + local_phase_order_values[q]) * 0.5 * fe_values.JxW(q);
              volume_1 +=
                (1 - local_phase_order_values[q]) * 0.5 * fe_values.JxW(q);
            }
        }
    }

  min_phase_value = Utilities::MPI::min(min_phase_value, mpi_communicator);
  max_phase_value = Utilities::MPI::max(max_phase_value, mpi_communicator);

  integral = Utilities::MPI::sum(integral, mpi_communicator);
  volume_0 = Utilities::MPI::sum(volume_0, mpi_communicator);
  volume_1 = Utilities::MPI::sum(volume_1, mpi_communicator);

  double global_volume = GridTools::volume(*triangulation, *mapping);
  double phase_average = integral / global_volume;


  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      announce_string(this->pcout, "Phase statistics");
      this->pcout << "Min: " << min_phase_value << std::endl;
      this->pcout << "Max: " << max_phase_value << std::endl;
      this->pcout << "Average: " << phase_average << std::endl;
      this->pcout << "Integral: " << integral << std::endl;
      this->pcout << "Volume phase 0: " << volume_0 << std::endl;
      this->pcout << "Volume phase 1: " << volume_1 << std::endl;
    }

  statistics_table.add_value("time", simulation_control->get_current_time());
  statistics_table.set_scientific("time", true);
  statistics_table.add_value("min", min_phase_value);
  statistics_table.set_scientific("min", true);
  statistics_table.add_value("max", max_phase_value);
  statistics_table.set_scientific("max", true);
  statistics_table.add_value("average", phase_average);
  statistics_table.set_scientific("average", true);
  statistics_table.add_value("integral", integral);
  statistics_table.set_scientific("integral", true);
  statistics_table.add_value("volume_0", volume_0);
  statistics_table.set_scientific("volume_0", true);
  statistics_table.add_value("volume_1", volume_1);
  statistics_table.set_scientific("volume_1", true);
}

template <int dim>
void
CahnHilliard<dim>::write_phase_statistics()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.phase_output_name + ".dat";
      std::ofstream output(filename.c_str());

      statistics_table.write_text(output);
    }
}

template <int dim>
void
CahnHilliard<dim>::calculate_phase_energy()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Scalar phase_order(0);

  const unsigned int          n_q_points = cell_quadrature->size();
  std::vector<double>         local_phase_order_values(n_q_points);
  std::vector<Tensor<1, dim>> local_phase_order_gradients(n_q_points);

  double bulk_energy(0.);
  double interface_energy(0.);
  double total_energy(0.);

  double epsilon =
    (this->simulation_parameters.multiphysics.cahn_hilliard_parameters
       .epsilon_set_method == Parameters::EpsilonSetMethod::manual) ?
      this->simulation_parameters.multiphysics.cahn_hilliard_parameters
        .epsilon :
      GridTools::minimal_cell_diameter(*triangulation);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[phase_order].get_function_values(present_solution,
                                                     local_phase_order_values);
          fe_values[phase_order].get_function_gradients(
            present_solution, local_phase_order_gradients);
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              bulk_energy += (1 - local_phase_order_values[q] *
                                    local_phase_order_values[q]) *
                             (1 - local_phase_order_values[q] *
                                    local_phase_order_values[q]) *
                             fe_values.JxW(q);
              interface_energy += epsilon * epsilon * 0.5 *
                                  (local_phase_order_gradients[q] *
                                   local_phase_order_gradients[q]) *
                                  fe_values.JxW(q);
            }
        }
    }

  bulk_energy      = Utilities::MPI::sum(bulk_energy, mpi_communicator);
  interface_energy = Utilities::MPI::sum(interface_energy, mpi_communicator);
  total_energy     = bulk_energy + interface_energy;

  phase_energy_table.add_value("time", simulation_control->get_current_time());
  phase_energy_table.set_scientific("time", true);
  phase_energy_table.add_value("bulk_energy", bulk_energy);
  phase_energy_table.set_scientific("bulk_energy", true);
  phase_energy_table.add_value("interface_energy", interface_energy);
  phase_energy_table.set_scientific("interface_energy", true);
  phase_energy_table.add_value("total_energy", total_energy);
  phase_energy_table.set_scientific("total_energy", true);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      announce_string(this->pcout, "Phase energy");
      this->pcout << "Bulk energy: "
                  << std::setprecision(simulation_control->get_log_precision())
                  << bulk_energy << std::endl;
      this->pcout << "Interface energy: "
                  << std::setprecision(simulation_control->get_log_precision())
                  << interface_energy << std::endl;
      this->pcout << "Total energy: "
                  << std::setprecision(simulation_control->get_log_precision())
                  << total_energy << std::endl;
    }
}

template <int dim>
void
CahnHilliard<dim>::write_phase_energy()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.phase_energy_output_name + ".dat";
      std::ofstream output(filename.c_str());
      phase_energy_table.set_precision("time", 12);
      phase_energy_table.set_precision("bulk_energy", 12);
      phase_energy_table.set_precision("interface_energy", 12);
      phase_energy_table.set_precision("total_energy", 12);

      phase_energy_table.write_text(output);
    }
}

template <int dim>
void
CahnHilliard<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      error_table.omit_column_from_convergence_rate_evaluation("cells");

      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_phase_order", true);
      error_table.set_precision("error_phase_order",
                                simulation_control->get_log_precision());
      error_table.set_scientific("error_chemical_potential", true);
      error_table.set_precision("error_chemical_potential",
                                simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
CahnHilliard<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
CahnHilliard<dim>::modify_solution()
{
  // Apply filter to phase order parameter
  apply_phase_filter();
}

template <int dim>
void
CahnHilliard<dim>::postprocess(bool first_iteration)
{
  auto         mpi_communicator = this->triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      double phase_order_error = calculate_L2_error().first;
      double potential_error   = calculate_L2_error().second;

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());
      error_table.add_value("error_phase_order", phase_order_error);
      error_table.add_value("error_chemical_potential", potential_error);

      error_table.set_scientific("error_phase_order", true);
      error_table.set_scientific("error_chemical_potential", true);

      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.analytical_solution->get_filename() +
        "_cahn_hilliard.dat";
      std::ofstream output(filename.c_str());
      error_table.write_text(output);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error phase order: " << phase_order_error
                      << std::endl;
          this->pcout << "L2 error chemical potential: " << potential_error
                      << std::endl;
        }
    }

  if (simulation_parameters.multiphysics.cahn_hilliard_parameters
        .epsilon_verbosity == Parameters::EpsilonVerbosity::verbose)
    {
      double epsilon = GridTools::minimal_cell_diameter(*triangulation);
      announce_string(this->pcout, "Epsilon value");
      this->pcout << "Epsilon value: " << epsilon << std::endl;
    }

  if (simulation_parameters.post_processing.calculate_phase_statistics)
    {
      calculate_phase_statistics();
      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_phase_statistics();
    }

  if (this->simulation_parameters.post_processing.calculate_phase_energy)
    {
      calculate_phase_energy();
      // Output phase energies to a text file from processor 0
      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        {
          this->write_phase_energy();
        }
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Cahn-Hilliard");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }

  if (this->simulation_parameters.post_processing.calculate_barycenter)
    {
      // Calculate volume and mass (this->mass_monitored)
      std::pair<Tensor<1, dim>, Tensor<1, dim>> position_and_velocity;

      if (multiphysics->fluid_dynamics_is_block())
        {
          // Check if the post processed variable needs to be calculated with
          // the average velocity profile or the fluid solution.
          if (this->simulation_parameters.initial_condition->type ==
                Parameters::FluidDynamicsInitialConditionType::
                  average_velocity_profile &&
              !this->simulation_parameters.multiphysics.fluid_dynamics &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing
                  .initial_time_for_average_velocities)
            {
              position_and_velocity = calculate_barycenter(
                this->present_solution,
                *multiphysics->get_block_time_average_solution(
                  PhysicsID::fluid_dynamics));
            }
          else
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_block_solution(
                                       PhysicsID::fluid_dynamics));
            }
        }
      else
        {
          // Check if the post processed variable needs to be calculated with
          // the average velocity profile or the fluid solution.
          if (this->simulation_parameters.initial_condition->type ==
                Parameters::FluidDynamicsInitialConditionType::
                  average_velocity_profile &&
              !this->simulation_parameters.multiphysics.fluid_dynamics &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing
                  .initial_time_for_average_velocities)
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_time_average_solution(
                                       PhysicsID::fluid_dynamics));
            }
          else
            {
              position_and_velocity =
                calculate_barycenter(this->present_solution,
                                     *multiphysics->get_solution(
                                       PhysicsID::fluid_dynamics));
            }
        }
      if (this_mpi_process == 0)
        {
          if (simulation_parameters.post_processing.verbosity ==
              Parameters::Verbosity::verbose)
            {
              std::cout << std::endl;
              std::string independent_column_names = "time";

              std::vector<std::string> dependent_column_names;
              dependent_column_names.emplace_back("x_cahn_hilliard");
              dependent_column_names.emplace_back("y_cahn_hilliard");
              if constexpr (dim == 3)
                dependent_column_names.emplace_back("z_cahn_hilliard");
              dependent_column_names.emplace_back("vx_cahn_hilliard");
              dependent_column_names.emplace_back("vy_cahn_hilliard");
              if constexpr (dim == 3)
                dependent_column_names.emplace_back("vz_cahn_hilliard");

              std::vector<Tensor<1, dim>> position_vector{
                position_and_velocity.first};
              std::vector<Tensor<1, dim>> velocity_vector{
                position_and_velocity.second};

              std::vector<std::vector<Tensor<1, dim>>>
                position_and_velocity_vectors{position_vector, velocity_vector};

              std::vector<double> time = {
                this->simulation_control->get_current_time()};

              TableHandler table = make_table_scalars_tensors(
                time,
                independent_column_names,
                position_and_velocity_vectors,
                dependent_column_names,
                this->simulation_parameters.simulation_control.log_precision,
                true);

              announce_string(this->pcout, "Cahn-Hilliard Barycenter");

              table.write_text(std::cout);
            }

          this->barycenter_table.add_value(
            "time", simulation_control->get_current_time());

          this->barycenter_table.add_value("x_cahn_hilliard",
                                           position_and_velocity.first[0]);
          this->barycenter_table.add_value("y_cahn_hilliard",
                                           position_and_velocity.first[1]);
          if constexpr (dim == 3)
            this->barycenter_table.add_value("z_cahn_hilliard",
                                             position_and_velocity.first[2]);

          this->barycenter_table.add_value("vx_cahn_hilliard",
                                           position_and_velocity.second[0]);
          this->barycenter_table.add_value("vy_cahn_hilliard",
                                           position_and_velocity.second[1]);
          if constexpr (dim == 3)
            this->barycenter_table.add_value("vz_cahn_hilliard",
                                             position_and_velocity.second[2]);


          if (this->simulation_control->get_step_number() %
                this->simulation_parameters.post_processing.output_frequency ==
              0)
            {
              // Save table to .dat
              std::string filename =
                this->simulation_parameters.simulation_control.output_folder +
                this->simulation_parameters.post_processing
                  .barycenter_output_name +
                ".dat";
              std::ofstream output(filename.c_str());
              this->barycenter_table.write_text(output);
              output.close();
            }
        }
    }
}


template <int dim>
void
CahnHilliard<dim>::pre_mesh_adaptation()
{
  solution_transfer->prepare_for_coarsening_and_refinement(present_solution);

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }
}

template <int dim>
void
CahnHilliard<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer
  GlobalVectorType tmp(locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer->interpolate(tmp);

  // Distribute constraints
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  // Transfer previous solutions
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      GlobalVectorType tmp_previous_solution(locally_owned_dofs,
                                             mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }
}

template <int dim>
void
CahnHilliard<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  const FEValuesExtractors::Scalar phase_order(0);
  const FEValuesExtractors::Scalar chemical_potential(1);

  if (ivar.first == Variable::phase_cahn_hilliard)
    {
      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(phase_order));
    }
  else if (ivar.first == Variable::chemical_potential_cahn_hilliard)
    {
      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(chemical_potential));
    }
}


template <int dim>
void
CahnHilliard<dim>::write_checkpoint()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(dof_handler);

  sol_set_transfer.emplace_back(&present_solution);
  for (const auto &previous_solution : previous_solutions)
    {
      sol_set_transfer.emplace_back(&previous_solution);
    }
  solution_transfer->prepare_for_serialization(sol_set_transfer);

  // Serialize tables
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    serialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_CH" + suffix,
      mpi_communicator);
  if (this->simulation_parameters.post_processing.calculate_phase_statistics)
    serialize_table(
      this->statistics_table,
      prefix + this->simulation_parameters.post_processing.phase_output_name +
        suffix,
      mpi_communicator);
}

template <int dim>
void
CahnHilliard<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading Cahn Hilliard checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions.size());
  GlobalVectorType distributed_system(locally_owned_dofs, mpi_communicator);
  input_vectors[0] = &distributed_system;


  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(locally_owned_dofs, mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  solution_transfer->deserialize(input_vectors);

  present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Deserialize tables
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    deserialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_CH" + suffix,
      mpi_communicator);
  if (this->simulation_parameters.post_processing.calculate_phase_statistics)
    deserialize_table(
      this->statistics_table,
      prefix + this->simulation_parameters.post_processing.phase_output_name +
        suffix,
      mpi_communicator);
}


template <int dim>
void
CahnHilliard<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  dof_handler.distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  auto mpi_communicator = triangulation->get_mpi_communicator();


  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  present_solution.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);

  filtered_solution.reinit(this->locally_owned_dofs,
                           this->locally_relevant_dofs,
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

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_cahn_hilliard.type)
      {
        ComponentMask mask(2, true);
        mask.set(1, false);

        // Dirichlet condition: imposed phase_order at i_bc
        // To impose the boundary condition only on the phase order, a component
        // mask is used at the end of the interpolate_boundary_values function
        if (type == BoundaryConditions::BoundaryType::
                      cahn_hilliard_dirichlet_phase_order)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              CahnHilliardFunctionDefined<dim>(
                &this->simulation_parameters.boundary_conditions_cahn_hilliard
                   .bcFunctions.at(id)
                   ->phi),
              nonzero_constraints,
              mask);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_cahn_hilliard
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_cahn_hilliard
                .periodic_direction.at(id),
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            zero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_cahn_hilliard.type)
      {
        if (type == BoundaryConditions::BoundaryType::
                      cahn_hilliard_dirichlet_phase_order)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              Functions::ZeroFunction<dim>(2),
              zero_constraints);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_cahn_hilliard
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_cahn_hilliard
                .periodic_direction.at(id),
              zero_constraints);
          }
      }
  }

  zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of Cahn-Hilliard degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the cahn_hilliard dof_handler and present solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::cahn_hilliard, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::cahn_hilliard, &this->present_solution);
  multiphysics->set_filtered_solution(PhysicsID::cahn_hilliard,
                                      &this->filtered_solution);
  multiphysics->set_previous_solutions(PhysicsID::cahn_hilliard,
                                       &this->previous_solutions);
}

template <int dim>
void
CahnHilliard<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions_cahn_hilliard
         .time_dependent)
    return;

  double time = this->simulation_control->get_current_time();
  nonzero_constraints.clear();
  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          nonzero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_cahn_hilliard.type)
    {
      ComponentMask mask(2, true);
      mask.set(1, false);

      // Dirichlet condition: imposed phase_order at i_bc
      // To impose the boundary condition only on the phase order, a component
      // mask is used at the end of the interpolate_boundary_values function
      if (type ==
          BoundaryConditions::BoundaryType::cahn_hilliard_dirichlet_phase_order)
        {
          this->simulation_parameters.boundary_conditions_cahn_hilliard
            .bcFunctions.at(id)
            ->phi.set_time(time);

          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            id,
            CahnHilliardFunctionDefined<dim>(
              &this->simulation_parameters.boundary_conditions_cahn_hilliard
                 .bcFunctions.at(id)
                 ->phi),
            nonzero_constraints,
            mask);
        }
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions_cahn_hilliard
              .periodic_neighbor_id.at(id),
            this->simulation_parameters.boundary_conditions_cahn_hilliard
              .periodic_direction.at(id),
            nonzero_constraints);
        }
    }

  nonzero_constraints.close();
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
CahnHilliard<dim>::set_initial_conditions()
{
  const FEValuesExtractors::Scalar phase_order(0);
  const FEValuesExtractors::Scalar potential(1);

  VectorTools::interpolate(
    *mapping,
    dof_handler,
    simulation_parameters.initial_condition->cahn_hilliard,
    newton_update,
    fe->component_mask(phase_order));

  // Set the initial chemical potential to 0. (May be discussed or modified
  // later)
  VectorTools::interpolate(
    *mapping,
    dof_handler,
    simulation_parameters.initial_condition->cahn_hilliard,
    newton_update,
    fe->component_mask(potential));

  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  apply_phase_filter();
  percolate_time_vectors();
}

template <int dim>
void
CahnHilliard<dim>::solve_linear_system(const bool initial_step,
                                       const bool /*renewed_matrix*/)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = triangulation->get_mpi_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const unsigned int ilu_fill = static_cast<unsigned int>(
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .ilu_precond_fill);
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(simulation_parameters.linear_solver
                                 .at(PhysicsID::cahn_hilliard)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
      .max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::cahn_hilliard)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template <int dim>
template <typename VectorType>
std::pair<Tensor<1, dim>, Tensor<1, dim>>
CahnHilliard<dim>::calculate_barycenter(const GlobalVectorType &solution,
                                        const VectorType       &solution_fd)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values_cahn_hilliard(*this->mapping,
                                        *this->fe,
                                        *this->cell_quadrature,
                                        update_values | update_gradients |
                                          update_quadrature_points |
                                          update_JxW_values);
  std::shared_ptr<CahnHilliardFilterBase> filter =
    CahnHilliardFilterBase::model_cast(
      this->simulation_parameters.multiphysics.cahn_hilliard_parameters);


  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  const unsigned int          n_q_points = this->cell_quadrature->size();
  std::vector<double>         phase_cahn_hilliard_values(n_q_points);
  std::vector<Tensor<1, dim>> phase_cahn_hilliard_gradients(n_q_points);
  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Point<dim>>     quadrature_locations(n_q_points);

  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar phase_order(0);

  Tensor<1, dim> barycenter_location;
  Tensor<1, dim> barycenter_velocity;
  double         volume = 0;


  std::map<field, std::vector<double>> fields;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_cahn_hilliard.reinit(cell);
          quadrature_locations =
            fe_values_cahn_hilliard.get_quadrature_points();
          fe_values_cahn_hilliard[phase_order].get_function_values(
            solution, phase_cahn_hilliard_values);
          fe_values_cahn_hilliard[phase_order].get_function_gradients(
            solution, phase_cahn_hilliard_gradients);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_fd);

          fe_values_fd.reinit(cell_fd);
          fe_values_fd[velocity].get_function_values(solution_fd,
                                                     velocity_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              const double JxW          = fe_values_cahn_hilliard.JxW(q);
              const double phase_values = phase_cahn_hilliard_values[q];



              volume += (1 - phase_values) * 0.5 * JxW;
              barycenter_location +=
                (1 - phase_values) * 0.5 * quadrature_locations[q] * JxW;
              barycenter_velocity +=
                (1 - phase_values) * 0.5 * velocity_values[q] * JxW;
            }
        }
    }

  volume = Utilities::MPI::sum(volume, mpi_communicator);
  barycenter_location =
    Utilities::MPI::sum(barycenter_location, mpi_communicator) / volume;
  barycenter_velocity =
    Utilities::MPI::sum(barycenter_velocity, mpi_communicator) / volume;

  return std::pair<Tensor<1, dim>, Tensor<1, dim>>(barycenter_location,
                                                   barycenter_velocity);
}

template <int dim>
void
CahnHilliard<dim>::apply_phase_filter()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  GlobalVectorType filtered_solution_owned(this->locally_owned_dofs,
                                           mpi_communicator);
  filtered_solution_owned = this->present_solution;

  filtered_solution.reinit(this->present_solution);

  // std::unordered_map<unsigned int, bool> filtered_cell_list;
  std::unordered_set<unsigned int> filtered_cell_list;

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Create filter object
  filter = CahnHilliardFilterBase::model_cast(
    this->simulation_parameters.multiphysics.cahn_hilliard_parameters);

  // Apply filter to solution
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int p = 0; p < local_dof_indices.size(); ++p)
            {
              if (this->locally_owned_dofs.is_element(local_dof_indices[p]))
                {
                  //  Allows to obtain the component corresponding to the degree
                  //  of freedom
                  auto component_index = fe->system_to_component_index(p).first;

                  // Filter only the phase field
                  if (component_index == 0)
                    {
                      auto iterator =
                        filtered_cell_list.find(local_dof_indices[p]);
                      if (iterator == filtered_cell_list.end())
                        {
                          filtered_cell_list.insert(local_dof_indices[p]);
                          filtered_solution_owned[local_dof_indices[p]] =
                            filter->filter_phase(
                              filtered_solution_owned[local_dof_indices[p]]);
                        }
                    }
                }
            }
        }
    }

  filtered_solution = filtered_solution_owned;

  if (this->simulation_parameters.multiphysics.cahn_hilliard_parameters
        .cahn_hilliard_phase_filter.verbosity == Parameters::Verbosity::verbose)
    {
      this->pcout << "Filtered phase values: " << std::endl;
      for (const double filtered_phase : filtered_solution)
        {
          this->pcout << filtered_phase << std::endl;
        }
    }
}

template <int dim>
void
CahnHilliard<dim>::output_newton_update_norms(
  const unsigned int display_precision)
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValuesExtractors::Scalar phase_order(0);
  FEValuesExtractors::Scalar chemical_potential(1);

  ComponentMask phase_order_mask = fe->component_mask(phase_order);
  ComponentMask chemical_potential_mask =
    fe->component_mask(chemical_potential);

  const std::vector<IndexSet> index_set_phase_order =
    DoFTools::locally_owned_dofs_per_component(dof_handler, phase_order_mask);
  const std::vector<IndexSet> index_set_chemical_potential =
    DoFTools::locally_owned_dofs_per_component(dof_handler,
                                               chemical_potential_mask);

  double local_sum = 0.0;
  double local_max = std::numeric_limits<double>::lowest();


  for (const auto &j : index_set_phase_order[0])
    {
      double dof_newton_update = newton_update[j];
      local_sum += dof_newton_update * dof_newton_update;
      local_max = std::max(local_max, std::abs(dof_newton_update));
    }


  double global_phase_order_l2_norm =
    std::sqrt(Utilities::MPI::sum(local_sum, mpi_communicator));
  double global_phase_order_linfty_norm =
    Utilities::MPI::max(local_max, mpi_communicator);

  local_sum = 0.0;
  local_max = std::numeric_limits<double>::lowest();

  for (const auto &j : index_set_chemical_potential[0])
    {
      double dof_newton_update = newton_update[j];
      local_sum += dof_newton_update * dof_newton_update;
      local_max = std::max(local_max, std::abs(dof_newton_update));
    }

  double global_chemical_potential_l2_norm =
    std::sqrt(Utilities::MPI::sum(local_sum, mpi_communicator));
  double global_chemical_potential_linfty_norm =
    Utilities::MPI::max(local_max, mpi_communicator);

  this->pcout << std::setprecision(display_precision)
              << "\n\t||dphi||_L2 = " << std::setw(6)
              << global_phase_order_l2_norm << std::setw(6)
              << "\t||dphi||_Linfty = " << std::setprecision(display_precision)
              << global_phase_order_linfty_norm << std::endl;
  this->pcout << std::setprecision(display_precision)
              << "\t||deta||_L2 = " << std::setw(6)
              << global_chemical_potential_l2_norm << std::setw(6)
              << "\t||deta||_Linfty = " << std::setprecision(display_precision)
              << global_chemical_potential_linfty_norm << std::endl;
}

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
CahnHilliard<2>::calculate_barycenter<GlobalVectorType>(
  const GlobalVectorType &solution,
  const GlobalVectorType &current_solution_fd);

template std::pair<Tensor<1, 3>, Tensor<1, 3>>
CahnHilliard<3>::calculate_barycenter<GlobalVectorType>(
  const GlobalVectorType &solution,
  const GlobalVectorType &current_solution_fd);

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
CahnHilliard<2>::calculate_barycenter<GlobalBlockVectorType>(
  const GlobalVectorType      &solution,
  const GlobalBlockVectorType &current_solution_fd);


template std::pair<Tensor<1, 3>, Tensor<1, 3>>
CahnHilliard<3>::calculate_barycenter<GlobalBlockVectorType>(
  const GlobalVectorType      &solution,
  const GlobalBlockVectorType &current_solution_fd);


template class CahnHilliard<2>;
template class CahnHilliard<3>;
