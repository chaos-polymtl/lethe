#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/vof.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_filter.h>
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
    this->simulation_control,
    this->simulation_parameters.fem_parameters,
    this->simulation_parameters.multiphysics.vof_parameters));
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
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics));
        }
    }
  else
    {
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics));
        }
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
  // TimerOutput::Scope t(this->computing_timer, "Assemble VOF RHS");

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
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_block_previous_solutions(
              PhysicsID::fluid_dynamics));
        }
    }
  else
    {
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            *multiphysics->get_previous_solutions(PhysicsID::fluid_dynamics));
        }
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

  // Filter phase fraction
  data_out.add_data_vector(this->dof_handler,
                           this->filtered_solution,
                           "filtered_phase");

  auto vof_parameters = this->simulation_parameters.multiphysics.vof_parameters;

  if (vof_parameters.peeling_wetting.enable_peeling ||
      vof_parameters.peeling_wetting.enable_wetting)
    {
      // Peeling/wetting output
      data_out.add_data_vector(this->dof_handler, marker_pw, "marker_pw");
    }

  if (vof_parameters.surface_tension_force.enable &&
      vof_parameters.surface_tension_force.output_vof_auxiliary_fields)
    {
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        projected_phase_fraction_gradient_component_interpretation(
          dim, DataComponentInterpretation::component_is_scalar);
      for (unsigned int i = 0; i < dim; ++i)
        projected_phase_fraction_gradient_component_interpretation[i] =
          DataComponentInterpretation::component_is_part_of_vector;

      std::vector<std::string> solution_names(dim, "phase_fraction_gradient");

      data_out.add_data_vector(
        projected_phase_fraction_gradient_dof_handler,
        present_projected_phase_fraction_gradient_solution,
        solution_names,
        projected_phase_fraction_gradient_component_interpretation);

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
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  const unsigned int  n_q_points = this->cell_quadrature->size();
  std::vector<double> phase_exact_solution(n_q_points);
  std::vector<double> phase_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->phase;
  exact_solution.set_time(this->simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          fe_values_vof.get_function_values(this->present_solution,
                                            phase_values);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values_vof.get_quadrature_points(),
                                    phase_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = phase_values[q];
              double exact = phase_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values_vof.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
template <typename VectorType>
std::pair<Tensor<1, dim>, Tensor<1, dim>>
VolumeOfFluid<dim>::calculate_barycenter(
  const TrilinosWrappers::MPI::Vector &solution,
  const VectorType &                   solution_fd)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);

  std::shared_ptr<VolumeOfFluidFilterBase> filter =
    VolumeOfFluidFilterBase::model_cast(
      this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  const unsigned int          n_q_points = this->cell_quadrature->size();
  std::vector<double>         phase_values(n_q_points);
  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Point<dim>>     quadrature_locations(n_q_points);

  const FEValuesExtractors::Vector velocity(0);

  Tensor<1, dim> barycenter_location;
  Tensor<1, dim> barycenter_velocity;
  double         volume = 0;


  std::map<field, std::vector<double>> fields;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          quadrature_locations = fe_values_vof.get_quadrature_points();
          fe_values_vof.get_function_values(solution, phase_values);

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
              const double JxW = fe_values_vof.JxW(q);
              const double filtered_phase_value =
                filter->filter_phase(phase_values[q]);

              volume += (filtered_phase_value)*JxW;
              barycenter_location +=
                (filtered_phase_value)*quadrature_locations[q] * JxW;
              barycenter_velocity +=
                (filtered_phase_value)*velocity_values[q] * JxW;
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

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
VolumeOfFluid<2>::calculate_barycenter<TrilinosWrappers::MPI::Vector>(
  const TrilinosWrappers::MPI::Vector &solution,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);


template std::pair<Tensor<1, 3>, Tensor<1, 3>>
VolumeOfFluid<3>::calculate_barycenter<TrilinosWrappers::MPI::Vector>(
  const TrilinosWrappers::MPI::Vector &solution,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template std::pair<Tensor<1, 2>, Tensor<1, 2>>
VolumeOfFluid<2>::calculate_barycenter<TrilinosWrappers::MPI::BlockVector>(
  const TrilinosWrappers::MPI::Vector &     solution,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);


template std::pair<Tensor<1, 3>, Tensor<1, 3>>
VolumeOfFluid<3>::calculate_barycenter<TrilinosWrappers::MPI::BlockVector>(
  const TrilinosWrappers::MPI::Vector &     solution,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);


template <int dim>
void
VolumeOfFluid<dim>::calculate_volume_and_mass(
  const TrilinosWrappers::MPI::Vector &solution,
  const Parameters::FluidIndicator     monitored_fluid)
{
  const MPI_Comm mpi_communicator = this->triangulation->get_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

  const unsigned int  n_q_points = this->cell_quadrature->size();
  std::vector<double> phase_values(n_q_points);
  std::vector<double> density_0(n_q_points);
  std::vector<double> density_1(n_q_points);

  this->volume_monitored = 0.;
  this->mass_monitored   = 0.;

  // Physical properties
  const auto density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::map<field, std::vector<double>> fields;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          fe_values_vof.get_function_values(solution, phase_values);

          // Calculate physical properties for the cell
          density_models[0]->vector_value(fields, density_0);
          density_models[1]->vector_value(fields, density_1);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (monitored_fluid)
                {
                    case Parameters::FluidIndicator::fluid0: {
                      this->volume_monitored +=
                        fe_values_vof.JxW(q) * (1 - phase_values[q]);
                      this->mass_monitored += fe_values_vof.JxW(q) *
                                              (1 - phase_values[q]) *
                                              density_0[q];
                      break;
                    }
                    case Parameters::FluidIndicator::fluid1: {
                      this->volume_monitored +=
                        fe_values_vof.JxW(q) * phase_values[q];
                      this->mass_monitored +=
                        fe_values_vof.JxW(q) * phase_values[q] * density_1[q];
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
                }
            }
        }
    }

  this->volume_monitored =
    Utilities::MPI::sum(this->volume_monitored, mpi_communicator);
  this->mass_monitored =
    Utilities::MPI::sum(this->mass_monitored, mpi_communicator);
}

template <int dim>
template <typename VectorType>
double
VolumeOfFluid<dim>::find_monitored_fluid_avg_pressure(
  const TrilinosWrappers::MPI::Vector &solution,
  const VectorType &                   current_solution_fd,
  const Parameters::FluidIndicator     monitored_fluid)
{
  QGauss<dim>    quadrature_formula(this->cell_quadrature->size());
  const MPI_Comm mpi_communicator = this->triangulation->get_communicator();

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values);

  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  FEValues<dim> fe_values_fd(*this->mapping,
                             dof_handler_fd->get_fe(),
                             quadrature_formula,
                             update_values);

  const unsigned int  n_q_points = quadrature_formula.size();
  std::vector<double> local_phase_values(n_q_points);

  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              local_pressure_values(n_q_points);

  double pressure_monitored_avg = 0.;
  int    nb_values(0);

  for (const auto &cell_vof : this->dof_handler.active_cell_iterators())
    {
      if (cell_vof->is_locally_owned())
        {
          fe_values_vof.reinit(cell_vof);
          fe_values_vof.get_function_values(solution, local_phase_values);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell_vof->level(),
            cell_vof->index(),
            dof_handler_fd);

          fe_values_fd.reinit(cell_fd);
          fe_values_fd[pressure].get_function_values(current_solution_fd,
                                                     local_pressure_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Gather minimum and maximum pressure on the phase of interest
              if ((monitored_fluid == Parameters::FluidIndicator::fluid0 &&
                   local_phase_values[q] < 0.5) ||
                  (monitored_fluid == Parameters::FluidIndicator::fluid1 &&
                   local_phase_values[q] > 0.5))
                {
                  pressure_monitored_avg += local_pressure_values[q];
                  nb_values++;
                }
            }
        }
    }

  pressure_monitored_avg =
    Utilities::MPI::sum(pressure_monitored_avg, mpi_communicator) /
    static_cast<double>(Utilities::MPI::sum(nb_values, mpi_communicator));

  return pressure_monitored_avg;
}

template double
VolumeOfFluid<2>::find_monitored_fluid_avg_pressure<
  TrilinosWrappers::MPI::Vector>(
  const TrilinosWrappers::MPI::Vector &solution,
  const TrilinosWrappers::MPI::Vector &current_solution_fd,
  const Parameters::FluidIndicator     monitored_fluid);

template double
VolumeOfFluid<3>::find_monitored_fluid_avg_pressure<
  TrilinosWrappers::MPI::Vector>(
  const TrilinosWrappers::MPI::Vector &solution,
  const TrilinosWrappers::MPI::Vector &current_solution_fd,
  const Parameters::FluidIndicator     monitored_fluid);

template double
VolumeOfFluid<2>::find_monitored_fluid_avg_pressure<
  TrilinosWrappers::MPI::BlockVector>(
  const TrilinosWrappers::MPI::Vector &     solution,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd,
  const Parameters::FluidIndicator          monitored_fluid);

template double
VolumeOfFluid<3>::find_monitored_fluid_avg_pressure<
  TrilinosWrappers::MPI::BlockVector>(
  const TrilinosWrappers::MPI::Vector &     solution,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd,
  const Parameters::FluidIndicator          monitored_fluid);

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
      if (this->simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        this->error_table.omit_column_from_convergence_rate_evaluation("cells");
      else
        this->error_table.omit_column_from_convergence_rate_evaluation("time");

      this->error_table.set_scientific("error_phase", true);
      this->error_table.set_precision(
        "error_phase", this->simulation_control->get_log_precision());
      this->error_table.write_text(std::cout);
    }
  if (this_mpi_process == 0 &&
      this->simulation_parameters.multiphysics.vof_parameters.conservation
          .verbosity == Parameters::Verbosity::extra_verbose)
    {
      this->table_monitoring_vof.write_text(std::cout);
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
VolumeOfFluid<dim>::postprocess(bool first_iteration)
{
  auto         mpi_communicator = this->triangulation->get_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

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

  if (this->simulation_parameters.multiphysics.vof_parameters.conservation
        .monitoring)
    {
      // Calculate volume and mass (this->mass_monitored)
      calculate_volume_and_mass(this->present_solution,
                                simulation_parameters.multiphysics
                                  .vof_parameters.conservation.monitored_fluid);

      if (first_iteration)
        this->mass_first_iteration = this->mass_monitored;



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

          std::string fluid_id("");

          if (this->simulation_parameters.multiphysics.vof_parameters
                .conservation.monitored_fluid ==
              Parameters::FluidIndicator::fluid1)
            {
              fluid_id = "fluid_1";
            }
          else if (this->simulation_parameters.multiphysics.vof_parameters
                     .conservation.monitored_fluid ==
                   Parameters::FluidIndicator::fluid0)
            {
              fluid_id = "fluid_0";
            }

          if (dim == 2)
            {
              // Add surface column
              this->table_monitoring_vof.add_value("surface_" + fluid_id,
                                                   this->volume_monitored);
              this->table_monitoring_vof.set_scientific("surface_" + fluid_id,
                                                        true);

              // Add mass per length column
              this->table_monitoring_vof.add_value("mass_per_length_" +
                                                     fluid_id,
                                                   this->mass_monitored);
              this->table_monitoring_vof.set_scientific("mass_per_length_" +
                                                          fluid_id,
                                                        true);
            }
          else if (dim == 3)
            {
              // Add volume column
              this->table_monitoring_vof.add_value("volume_" + fluid_id,
                                                   this->volume_monitored);
              this->table_monitoring_vof.set_scientific("volume_" + fluid_id,
                                                        true);

              // Add mass column
              this->table_monitoring_vof.add_value("mass_" + fluid_id,
                                                   this->mass_monitored);
              this->table_monitoring_vof.set_scientific("mass_" + fluid_id,
                                                        true);
            }

          // Add sharpening threshold column
          this->table_monitoring_vof.add_value("sharpening_threshold",
                                               this->sharpening_threshold);

          if (this->simulation_control->get_step_number() %
                this->simulation_parameters.post_processing.output_frequency ==
              0)
            {
              // Save table to .dat
              std::string filename =
                this->simulation_parameters.simulation_control.output_folder +
                "VOF_monitoring_" + fluid_id + ".dat";
              std::ofstream output(filename.c_str());
              this->table_monitoring_vof.write_text(output);
            }
        }
    }

  if (this->simulation_parameters.post_processing.calculate_vof_barycenter)
    {
      // Calculate volume and mass (this->mass_monitored)
      std::pair<Tensor<1, dim>, Tensor<1, dim>> position_and_velocity;

      if (multiphysics->fluid_dynamics_is_block())
        {
          if (this->simulation_parameters.multiphysics
                .use_time_average_velocity_field &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing.initial_time)
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
          if (this->simulation_parameters.multiphysics
                .use_time_average_velocity_field &&
              simulation_control->get_current_time() >
                this->simulation_parameters.post_processing.initial_time)
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
              dependent_column_names.push_back("x_vof");
              dependent_column_names.push_back("y_vof");
              if (dim == 3)
                dependent_column_names.push_back("z_vof");
              dependent_column_names.push_back("vx_vof");
              dependent_column_names.push_back("vy_vof");
              if (dim == 3)
                dependent_column_names.push_back("vz_vof");

              std::vector<Tensor<1, dim>> position_and_velocity_vector;
              position_and_velocity_vector.push_back(
                position_and_velocity.first);
              position_and_velocity_vector.push_back(
                position_and_velocity.second);

              std::vector<double> time(
                this->simulation_control->get_current_time());

              TableHandler table = make_table_scalars_tensors(
                time,
                independent_column_names,
                position_and_velocity_vector,
                dependent_column_names,
                this->simulation_parameters.simulation_control.log_precision);

              std::cout << "+------------------------------------------+"
                        << std::endl;
              std::cout << "|  VOF Barycenter                          |"
                        << std::endl;
              std::cout << "+------------------------------------------+"
                        << std::endl;
              table.write_text(std::cout);
            }

          this->table_barycenter.add_value(
            "time", simulation_control->get_current_time());

          this->table_barycenter.add_value("x_vof",
                                           position_and_velocity.first[0]);
          this->table_barycenter.add_value("y_vof",
                                           position_and_velocity.first[1]);
          if constexpr (dim == 3)
            this->table_barycenter.add_value("z_vof",
                                             position_and_velocity.first[2]);

          this->table_barycenter.add_value("vx_vof",
                                           position_and_velocity.second[0]);
          this->table_barycenter.add_value("vy_vof",
                                           position_and_velocity.second[1]);
          if constexpr (dim == 3)
            this->table_barycenter.add_value("vz_vof",
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
              this->table_barycenter.write_text(output);
              output.close();
            }
        }
    }
}

template <int dim>
void
VolumeOfFluid<dim>::modify_solution()
{
  auto vof_parameters = this->simulation_parameters.multiphysics.vof_parameters;
  // Peeling/wetting
  if (vof_parameters.peeling_wetting.enable_peeling ||
      vof_parameters.peeling_wetting.enable_wetting)
    {
      handle_peeling_wetting();
    }
  // Interface sharpening
  if (vof_parameters.sharpening.enable)
    {
      // Interface sharpening is done at a constant frequency
      if (this->simulation_control->get_step_number() %
            this->simulation_parameters.multiphysics.vof_parameters.sharpening
              .frequency ==
          0)
        {
          handle_interface_sharpening();
        }
    }
  // Apply filter to phase fraction
  apply_phase_filter();
  if (vof_parameters.surface_tension_force.enable)
    {
      find_projected_phase_fraction_gradient();
      find_projected_interface_curvature();
    }
}

template <int dim>
void
VolumeOfFluid<dim>::handle_interface_sharpening()
{
  if (this->simulation_parameters.multiphysics.vof_parameters.sharpening
          .verbosity != Parameters::Verbosity::quiet ||
      this->simulation_parameters.multiphysics.vof_parameters.conservation
          .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "Sharpening interface at step "
                  << this->simulation_control->get_step_number() << std::endl;
    }
  if (this->simulation_parameters.multiphysics.vof_parameters.sharpening.type ==
      Parameters::SharpeningType::adaptative)
    {
      if (this->simulation_parameters.multiphysics.vof_parameters.conservation
            .verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << "   Adapting the sharpening threshold" << std::endl;
        }

      this->sharpening_threshold = find_sharpening_threshold();

      if (this->simulation_parameters.multiphysics.vof_parameters.conservation
            .verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << "   ... final sharpening is : "
                      << this->sharpening_threshold << std::endl;
        }
    }
  else
    {
      // Constant sharpening
      this->sharpening_threshold = this->simulation_parameters.multiphysics
                                     .vof_parameters.sharpening.threshold;
    }

  // Sharpen the interface of all solutions (present and previous)
  sharpen_interface(this->present_solution, this->sharpening_threshold, true);
}

template <int dim>
double
VolumeOfFluid<dim>::find_sharpening_threshold()
{
  // Sharpening threshold (st) search range extrema
  double st_min = 0.5 - this->simulation_parameters.multiphysics.vof_parameters
                          .sharpening.threshold_max_deviation;
  double st_max = 0.5 + this->simulation_parameters.multiphysics.vof_parameters
                          .sharpening.threshold_max_deviation;

  // Useful definitions for readability
  const double mass_deviation_tol = this->simulation_parameters.multiphysics
                                      .vof_parameters.conservation.tolerance *
                                    this->mass_first_iteration;
  const unsigned int max_iterations =
    this->simulation_parameters.multiphysics.vof_parameters.sharpening
      .max_iterations;

  const Parameters::FluidIndicator monitored_fluid =
    this->simulation_parameters.multiphysics.vof_parameters.conservation
      .monitored_fluid;

  unsigned int nb_search_ite = 0;
  // Local variable for the tested sharpening_threshold values

  double mass_deviation_min = calculate_mass_deviation(monitored_fluid, st_min);
  double mass_deviation_max = calculate_mass_deviation(monitored_fluid, st_max);
  double mass_deviation_avg = 0.;
  double st_avg             = 0.;


  // Bissection algorithm to calculate an interface sharpening value that would
  // ensure mass conservation of the monitored phase (do-while loop, see
  // condition below)
  do
    {
      nb_search_ite++;
      // Calculate middle point value
      st_avg = (st_min + st_max) / 2.;

      mass_deviation_avg = calculate_mass_deviation(monitored_fluid, st_avg);

      if (this->simulation_parameters.multiphysics.vof_parameters.conservation
            .verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout
            << "   ... step " << nb_search_ite
            << " of the search algorithm, min, avg, max mass deviation is : "
            << mass_deviation_min << " " << mass_deviation_avg << " "
            << mass_deviation_max << std::endl;
        }

      if (mass_deviation_avg * mass_deviation_min < 0)
        {
          st_max             = st_avg;
          mass_deviation_max = mass_deviation_avg;
        }
      else
        {
          st_min             = st_avg;
          mass_deviation_min = mass_deviation_avg;
        }
  } while (std::abs(mass_deviation_avg) > mass_deviation_tol &&
           nb_search_ite < max_iterations);

  // Take minimum deviation in between the two endpoints of the last
  // interval searched, if out of the do-while loop because max_iterations is
  // reached
  if (std::abs(mass_deviation_avg) > mass_deviation_tol)
    {
      double mass_deviation_endpoint = 0.;
      double st_tested               = 0;
      if (st_min == st_avg)
        st_tested = st_max;
      else if (st_max == st_avg)
        st_tested = st_min;

      mass_deviation_endpoint =
        calculate_mass_deviation(monitored_fluid, st_tested);

      // Retake st_ave value if mass deviation is not lowered at endpoint
      // values
      if (std::abs(mass_deviation_endpoint) > std::abs(mass_deviation_avg))
        {
          st_tested = st_avg;
        }

      // Output message
      if (std::abs(mass_deviation_endpoint) > mass_deviation_tol)
        {
          this->pcout
            << "  WARNING: Maximum number of iterations (" << nb_search_ite
            << ") reached in the " << std::endl
            << "  adaptative sharpening threshold algorithm, remaining error"
            << std::endl
            << "  on mass conservation is: "
            << (this->mass_monitored - this->mass_first_iteration) /
                 this->mass_first_iteration
            << std::endl
            << "  Consider increasing the sharpening threshold range or the "
            << std::endl
            << "  number of iterations to reach the mass conservation tolerance."
            << std::endl;
        }
    }

  // Output message that mass conservation condition is reached
  if (this->simulation_parameters.multiphysics.vof_parameters.conservation
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "   ... search algorithm took : " << nb_search_ite
                  << " step(s) " << std::endl
                  << "   ... error on mass conservation reached: "
                  << (this->mass_monitored - this->mass_first_iteration) /
                       this->mass_first_iteration
                  << std::endl;
    }

  return st_avg;
}

template <int dim>
double
VolumeOfFluid<dim>::calculate_mass_deviation(
  const Parameters::FluidIndicator monitored_fluid,
  const double                     sharpening_threshold)
{
  // Copy present solution VOF
  evaluation_point = this->present_solution;

  // Sharpen interface using the tested threshold value
  sharpen_interface(evaluation_point, sharpening_threshold, false);

  // Calculate mass of the monitored phase
  calculate_volume_and_mass(evaluation_point, monitored_fluid);

  // Calculate relative mass deviation
  double mass_deviation = (this->mass_monitored - this->mass_first_iteration) /
                          this->mass_first_iteration;

  return mass_deviation;
}

template <int dim>
void
VolumeOfFluid<dim>::sharpen_interface(TrilinosWrappers::MPI::Vector &solution,
                                      const double sharpening_threshold,
                                      const bool   sharpen_previous_solutions)
{
  // Limit the phase fractions between 0 and 1
  update_solution_and_constraints(solution);
  if (sharpen_previous_solutions)
    {
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        update_solution_and_constraints(previous_solutions[p]);
    }

  // Assemble matrix and solve the system for interface sharpening
  assemble_L2_projection_interface_sharpening(solution, sharpening_threshold);
  solve_interface_sharpening(solution);

  if (sharpen_previous_solutions)
    {
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        {
          assemble_L2_projection_interface_sharpening(previous_solutions[p],
                                                      sharpening_threshold);
          solve_interface_sharpening(previous_solutions[p]);
        }
    }

  // Re limit the phase fractions between 0 and 1 after interface
  // sharpening
  update_solution_and_constraints(solution);
  if (sharpen_previous_solutions)
    {
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        update_solution_and_constraints(previous_solutions[p]);
    }
}

template <int dim>
void
VolumeOfFluid<dim>::smooth_phase_fraction()
{
  assemble_projection_phase_fraction(present_solution);
  solve_projection_phase_fraction(present_solution);
}

template <int dim>
void
VolumeOfFluid<dim>::find_projected_phase_fraction_gradient()
{
  assemble_projected_phase_fraction_gradient_matrix_and_rhs(filtered_solution);
  solve_projected_phase_fraction_gradient();
}

template <int dim>
void
VolumeOfFluid<dim>::find_projected_interface_curvature()
{
  assemble_curvature_matrix_and_rhs(
    present_projected_phase_fraction_gradient_solution);
  solve_curvature();
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_projection_phase_fraction(
  TrilinosWrappers::MPI::Vector &solution)
{
  // Get fe values of VOF phase fraction
  FEValues<dim> fe_values_phase_fraction(*this->mapping,
                                         *this->fe,
                                         *this->cell_quadrature,
                                         update_values | update_JxW_values |
                                           update_gradients);

  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_phase_fraction(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_phase_fraction(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double>         phi_phase_fraction(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_phase_fraction_gradient(dofs_per_cell);

  std::vector<double> phase_values(n_q_points);

  // Reinitialize system matrix and rhs for the projected phase fraction
  system_rhs_phase_fraction    = 0;
  system_matrix_phase_fraction = 0;

  double h;
  double cell_volume;

  const double phase_fraction_diffusion_factor =
    this->simulation_parameters.initial_condition
      ->projection_step_diffusion_factor;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_phase_fraction.reinit(cell);

          cell_volume = cell->measure();
          if (dim == 2)
            h = std::sqrt(4. * cell_volume / M_PI) / fe->degree;
          else if (dim == 3)
            h = pow(6 * cell_volume / M_PI, 1. / 3.) / fe->degree;

          local_matrix_phase_fraction = 0;
          local_rhs_phase_fraction    = 0;

          // Get phase fraction values
          fe_values_phase_fraction.get_function_values(solution, phase_values);

          double color_function = 0.0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              color_function +=
                phase_values[q] * fe_values_phase_fraction.JxW(q);
            }

          color_function /= cell_volume;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_phase_fraction[k] =
                    fe_values_phase_fraction.shape_value(k, q);
                  phi_phase_fraction_gradient[k] =
                    fe_values_phase_fraction.shape_grad(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_phase_fraction(i, j) +=
                        (phi_phase_fraction[j] * phi_phase_fraction[i] +
                         phase_fraction_diffusion_factor * h * h *
                           scalar_product(phi_phase_fraction_gradient[i],
                                          phi_phase_fraction_gradient[j])) *
                        fe_values_phase_fraction.JxW(q);
                    }

                  // rhs
                  local_rhs_phase_fraction(i) +=
                    phi_phase_fraction[i] * color_function *
                    fe_values_phase_fraction.JxW(q);
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
VolumeOfFluid<dim>::solve_projection_phase_fraction(
  TrilinosWrappers::MPI::Vector &solution)
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-13;

  TrilinosWrappers::MPI::Vector completely_distributed_phase_fraction_solution(
    this->locally_owned_dofs, triangulation->get_communicator());

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

  ilu_preconditioner->initialize(system_matrix_phase_fraction,
                                 preconditionerOptions);
  solver.solve(system_matrix_phase_fraction,
               completely_distributed_phase_fraction_solution,
               system_rhs_phase_fraction,
               *ilu_preconditioner);

  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver (phase fraction gradient) took : "
                  << solver_control.last_step() << " steps " << std::endl;
    }

  this->nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_projected_phase_fraction_gradient_matrix_and_rhs(
  TrilinosWrappers::MPI::Vector &solution)
{
  // Get fe values of VOF phase fraction and phase fraction gradient (pfg)
  FEValues<dim> fe_values_phase_fraction(*this->mapping,
                                         *this->fe,
                                         *this->cell_quadrature,
                                         update_values | update_gradients);

  FEValues<dim> fe_values_projected_phase_fraction_gradient(
    *this->mapping,
    *this->fe_projected_phase_fraction_gradient,
    *this->cell_quadrature,
    update_values | update_JxW_values | update_gradients);


  // const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int dofs_per_cell =
    this->fe_projected_phase_fraction_gradient->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_projected_phase_fraction_gradient(
    dofs_per_cell, dofs_per_cell);
  Vector<double> local_rhs_projected_phase_fraction_gradient(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<Tensor<1, dim>> phi_projected_phase_fraction_gradient(
    dofs_per_cell);
  std::vector<Tensor<2, dim>> phi_projected_phase_fraction_gradient_gradient(
    dofs_per_cell);

  const FEValuesExtractors::Vector pfg(0);

  std::vector<Tensor<1, dim>> phase_gradient_values(n_q_points);
  std::vector<Tensor<1, dim>> projected_phase_fraction_gradient_values(
    n_q_points);
  std::vector<Tensor<2, dim>> projected_phase_fraction_gradient_gradient_values(
    n_q_points);

  // Reinitialize system matrix and rhs for the pfg
  system_rhs_projected_phase_fraction_gradient    = 0;
  system_matrix_projected_phase_fraction_gradient = 0;

  double h;

  const double phase_fraction_gradient_diffusion_factor =
    this->simulation_parameters.multiphysics.vof_parameters
      .surface_tension_force.phase_fraction_gradient_diffusion_factor;

  for (const auto &projected_phase_fraction_gradient_cell :
       this->projected_phase_fraction_gradient_dof_handler
         .active_cell_iterators())
    {
      if (projected_phase_fraction_gradient_cell->is_locally_owned())
        {
          // Gather the active cell iterator related to the VOF phase fraction
          typename DoFHandler<dim>::active_cell_iterator cell(
            &(*this->triangulation),
            projected_phase_fraction_gradient_cell->level(),
            projected_phase_fraction_gradient_cell->index(),
            &this->dof_handler);

          fe_values_phase_fraction.reinit(cell);
          fe_values_projected_phase_fraction_gradient.reinit(
            projected_phase_fraction_gradient_cell);

          auto &fe_filtered_phase_fraction =
            fe_values_projected_phase_fraction_gradient.get_fe();

          if (dim == 2)
            h =
              std::sqrt(4. * projected_phase_fraction_gradient_cell->measure() /
                        M_PI) /
              fe_filtered_phase_fraction.degree;
          else if (dim == 3)
            h =
              pow(6 * projected_phase_fraction_gradient_cell->measure() / M_PI,
                  1. / 3.) /
              fe_filtered_phase_fraction.degree;


          local_matrix_projected_phase_fraction_gradient = 0;
          local_rhs_projected_phase_fraction_gradient    = 0;

          // Get phase fraction values, pfg values and gradients
          fe_values_phase_fraction.get_function_gradients(
            solution, phase_gradient_values);

          fe_values_projected_phase_fraction_gradient[pfg].get_function_values(
            present_projected_phase_fraction_gradient_solution,
            projected_phase_fraction_gradient_values);

          fe_values_projected_phase_fraction_gradient[pfg]
            .get_function_gradients(
              present_projected_phase_fraction_gradient_solution,
              projected_phase_fraction_gradient_gradient_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_projected_phase_fraction_gradient[k] =
                    fe_values_projected_phase_fraction_gradient[pfg].value(k,
                                                                           q);
                  phi_projected_phase_fraction_gradient_gradient[k] =
                    fe_values_projected_phase_fraction_gradient[pfg].gradient(
                      k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_projected_phase_fraction_gradient(i, j) +=
                        (phi_projected_phase_fraction_gradient[j] *
                           phi_projected_phase_fraction_gradient[i] +
                         h * h * phase_fraction_gradient_diffusion_factor *
                           scalar_product(
                             phi_projected_phase_fraction_gradient_gradient[i],
                             phi_projected_phase_fraction_gradient_gradient
                               [j])) *
                        fe_values_projected_phase_fraction_gradient.JxW(q);
                    }

                  // rhs
                  local_rhs_projected_phase_fraction_gradient(i) +=
                    phi_projected_phase_fraction_gradient[i] *
                    phase_gradient_values[q] *
                    fe_values_projected_phase_fraction_gradient.JxW(q);
                }
            }

          projected_phase_fraction_gradient_cell->get_dof_indices(
            local_dof_indices);
          projected_phase_fraction_gradient_constraints
            .distribute_local_to_global(
              local_matrix_projected_phase_fraction_gradient,
              local_rhs_projected_phase_fraction_gradient,
              local_dof_indices,
              system_matrix_projected_phase_fraction_gradient,
              system_rhs_projected_phase_fraction_gradient);
        }
    }
  system_matrix_projected_phase_fraction_gradient.compress(
    VectorOperation::add);
  system_rhs_projected_phase_fraction_gradient.compress(VectorOperation::add);
}


template <int dim>
void
VolumeOfFluid<dim>::solve_projected_phase_fraction_gradient()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-13;

  TrilinosWrappers::MPI::Vector
    completely_distributed_projected_phase_fraction_gradient_solution(
      this->locally_owned_dofs_projected_phase_fraction_gradient,
      triangulation->get_communicator());

  completely_distributed_projected_phase_fraction_gradient_solution =
    present_projected_phase_fraction_gradient_solution;

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

  ilu_preconditioner->initialize(
    system_matrix_projected_phase_fraction_gradient, preconditionerOptions);
  solver.solve(
    system_matrix_projected_phase_fraction_gradient,
    completely_distributed_projected_phase_fraction_gradient_solution,
    system_rhs_projected_phase_fraction_gradient,
    *ilu_preconditioner);

  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver (phase fraction gradient) took : "
                  << solver_control.last_step() << " steps " << std::endl;
    }

  projected_phase_fraction_gradient_constraints.distribute(
    completely_distributed_projected_phase_fraction_gradient_solution);
  present_projected_phase_fraction_gradient_solution =
    completely_distributed_projected_phase_fraction_gradient_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_curvature_matrix_and_rhs(
  TrilinosWrappers::MPI::Vector
    &present_projected_phase_fraction_gradient_solution)
{
  // Get fe values of phase fraction gradient (pfg) and curvature
  FEValues<dim> fe_values_curvature(*this->curvature_mapping,
                                    *this->fe_curvature,
                                    *this->cell_quadrature,
                                    update_values | update_JxW_values |
                                      update_gradients);

  FEValues<dim> fe_values_projected_phase_fraction_gradient(
    *this->mapping,
    *this->fe_projected_phase_fraction_gradient,
    *this->cell_quadrature,
    update_values);

  const unsigned int dofs_per_cell = this->fe_curvature->dofs_per_cell;

  const unsigned int n_q_points = this->cell_quadrature->size();
  FullMatrix<double> local_matrix_curvature(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs_curvature(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double>         phi_curvature(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_curvature_gradient(dofs_per_cell);

  projected_phase_fraction_gradient_values =
    std::vector<Tensor<1, dim>>(n_q_points);
  curvature_values = std::vector<double>(n_q_points);
  std::vector<Tensor<1, dim>>      curvature_gradient_values(n_q_points);
  const FEValuesExtractors::Vector pfg(0);

  // Reinitialize system matrix and rhs for the curvature
  system_rhs_curvature    = 0;
  system_matrix_curvature = 0;

  double curvature_diffusion_factor =
    simulation_parameters.multiphysics.vof_parameters.surface_tension_force
      .curvature_diffusion_factor;

  double h;

  for (const auto &curvature_cell :
       this->curvature_dof_handler.active_cell_iterators())
    {
      if (curvature_cell->is_locally_owned())
        {
          // Gather the active cell iterator related to the phase fraction
          // gradient (pfg)
          typename DoFHandler<dim>::active_cell_iterator
            projected_phase_fraction_gradient_cell(
              &(*this->triangulation),
              curvature_cell->level(),
              curvature_cell->index(),
              &this->projected_phase_fraction_gradient_dof_handler);

          fe_values_projected_phase_fraction_gradient.reinit(
            projected_phase_fraction_gradient_cell);
          fe_values_curvature.reinit(curvature_cell);

          local_matrix_curvature = 0;
          local_rhs_curvature    = 0;

          auto &fe_curvature = fe_values_curvature.get_fe();

          if (dim == 2)
            h = std::sqrt(4. * curvature_cell->measure() / M_PI) /
                fe_curvature.degree;
          else if (dim == 3)
            h = pow(6 * curvature_cell->measure() / M_PI, 1. / 3.) /
                fe_curvature.degree;

          // Get pfg values, curvature values and gradients
          fe_values_projected_phase_fraction_gradient[pfg].get_function_values(
            present_projected_phase_fraction_gradient_solution,
            projected_phase_fraction_gradient_values);

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

              // Calculate phase gradient norm and add a tolerance to it to
              // prevent illegal divisions
              const double projected_phase_fraction_gradient_norm =
                projected_phase_fraction_gradient_values[q].norm() + 1e-12;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix_curvature(i, j) +=
                        (phi_curvature[j] * phi_curvature[i] +
                         h * h * curvature_diffusion_factor *
                           scalar_product(phi_curvature_gradient[i],
                                          phi_curvature_gradient[j])) *
                        fe_values_curvature.JxW(q);
                    }
                  // rhs

                  local_rhs_curvature(i) +=
                    phi_curvature_gradient[i] *
                    (projected_phase_fraction_gradient_values[q] /
                     projected_phase_fraction_gradient_norm) *
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
  const double linear_solver_tolerance = 1e-13;

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

  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.verbosity != Parameters::Verbosity::quiet)
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
  this->solution_transfer->prepare_for_coarsening_and_refinement(
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
  this->solution_transfer->interpolate(tmp);

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

  // PFG and curvature
  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.enable)
    {
      find_projected_phase_fraction_gradient();
      find_projected_interface_curvature();
    }
}

template <int dim>
void
VolumeOfFluid<dim>::compute_kelly(
  const std::pair<const Parameters::MeshAdaptation::Variable,
                  Parameters::MultipleAdaptationParameters> &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Parameters::MeshAdaptation::Variable::phase)
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

  solution_transfer =
    std::make_shared<parallel::distributed::
                       SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>(
      dof_handler);

  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < this->previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&this->previous_solutions[i]);
    }
  this->solution_transfer->prepare_for_serialization(sol_set_transfer);
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

  this->solution_transfer->deserialize(input_vectors);

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

  if (this->simulation_parameters.multiphysics.vof_parameters
        .surface_tension_force.enable)
    {
      projected_phase_fraction_gradient_dof_handler.distribute_dofs(
        *fe_projected_phase_fraction_gradient);

      locally_owned_dofs_projected_phase_fraction_gradient =
        projected_phase_fraction_gradient_dof_handler.locally_owned_dofs();

      DoFTools::extract_locally_relevant_dofs(
        projected_phase_fraction_gradient_dof_handler,
        locally_relevant_dofs_projected_phase_fraction_gradient);

      projected_phase_fraction_gradient_constraints.clear();
      projected_phase_fraction_gradient_constraints.reinit(
        locally_relevant_dofs_projected_phase_fraction_gradient);
      DoFTools::make_hanging_node_constraints(
        projected_phase_fraction_gradient_dof_handler,
        projected_phase_fraction_gradient_constraints);
      projected_phase_fraction_gradient_constraints.close();

      nodal_projected_phase_fraction_gradient_relevant.reinit(
        locally_owned_dofs_projected_phase_fraction_gradient,
        locally_relevant_dofs_projected_phase_fraction_gradient,
        mpi_communicator);

      nodal_projected_phase_fraction_gradient_owned.reinit(
        locally_owned_dofs_projected_phase_fraction_gradient, mpi_communicator);



      DynamicSparsityPattern dsp_projected_phase_fraction_gradient(
        locally_relevant_dofs_projected_phase_fraction_gradient);
      DoFTools::make_sparsity_pattern(
        projected_phase_fraction_gradient_dof_handler,
        dsp_projected_phase_fraction_gradient,
        projected_phase_fraction_gradient_constraints,
        false);

      SparsityTools::distribute_sparsity_pattern(
        dsp_projected_phase_fraction_gradient,
        locally_owned_dofs_projected_phase_fraction_gradient,
        mpi_communicator,
        locally_relevant_dofs_projected_phase_fraction_gradient);


      // Initialization of phase fraction gradient matrix and rhs for the
      // calculation surface tension force
      system_matrix_projected_phase_fraction_gradient.reinit(
        locally_owned_dofs_projected_phase_fraction_gradient,
        locally_owned_dofs_projected_phase_fraction_gradient,
        dsp_projected_phase_fraction_gradient,
        mpi_communicator);

      system_rhs_projected_phase_fraction_gradient.reinit(
        locally_owned_dofs_projected_phase_fraction_gradient, mpi_communicator);

      present_projected_phase_fraction_gradient_solution.reinit(
        locally_owned_dofs_projected_phase_fraction_gradient,
        locally_relevant_dofs_projected_phase_fraction_gradient,
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

      nodal_curvature_owned.reinit(locally_owned_dofs_curvature,
                                   mpi_communicator);



      DynamicSparsityPattern dsp_curvature(locally_relevant_dofs_curvature);
      DoFTools::make_sparsity_pattern(curvature_dof_handler,
                                      dsp_curvature,
                                      curvature_constraints,
                                      false);

      SparsityTools::distribute_sparsity_pattern(
        dsp_curvature,
        locally_owned_dofs_curvature,
        mpi_communicator,
        locally_relevant_dofs_curvature);


      // Initialization of curvature matrix and rhs for the
      // calculation surface tension force
      system_matrix_curvature.reinit(locally_owned_dofs_curvature,
                                     locally_owned_dofs_curvature,
                                     dsp_curvature,
                                     mpi_communicator);

      system_rhs_curvature.reinit(locally_owned_dofs_curvature,
                                  mpi_communicator);

      present_curvature_solution.reinit(locally_owned_dofs_curvature,
                                        locally_relevant_dofs_curvature,
                                        mpi_communicator);
    }


  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                mpi_communicator);

  this->filtered_solution.reinit(this->locally_owned_dofs,
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

  define_non_zero_constraints();

  // Boundary conditions for Newton correction
  define_zero_constraints();

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
  // used in update_solution_and_constraints to limit the phase fraction
  // values between 0 and 1. Accoring to step-41, to limit the phase fractions
  // we compute the Lagrange multiplier as the residual of the original linear
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

  this->pcout << "   Number of VOF degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;

  // Provide the VOF dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::VOF, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::VOF, &this->present_solution);
  multiphysics->set_filtered_solution(PhysicsID::VOF, &this->filtered_solution);
  multiphysics->set_previous_solutions(PhysicsID::VOF,
                                       &this->previous_solutions);


  mass_matrix_phase_fraction.reinit(this->locally_owned_dofs,
                                    this->locally_owned_dofs,
                                    dsp,
                                    mpi_communicator);

  assemble_mass_matrix(mass_matrix_phase_fraction);
}

template <int dim>
void
VolumeOfFluid<dim>::define_zero_constraints()
{
  // Zero constraints
  this->zero_constraints.clear();
  this->zero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
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
    }
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_vof.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions_vof.type[i_bc] ==
          BoundaryConditions::BoundaryType::vof_dirichlet)
        {
          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            this->simulation_parameters.boundary_conditions_vof.id[i_bc],
            Functions::ZeroFunction<dim>(),
            this->zero_constraints);
        }
    }
  this->zero_constraints.close();
}

template <int dim>
void
VolumeOfFluid<dim>::define_non_zero_constraints()
{
  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(this->locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);

    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
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
    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions_vof.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions_vof.type[i_bc] ==
            BoundaryConditions::BoundaryType::vof_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_vof.id[i_bc],
              *this->simulation_parameters.boundary_conditions_vof
                 .phase_fraction[i_bc],
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();
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
  apply_phase_filter();

  if (simulation_parameters.initial_condition->enable_projection_step)
    smooth_phase_fraction();

  percolate_time_vectors();
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
  TrilinosWrappers::MPI::Vector &solution,
  const double                   sharpening_threshold)
{
  const double interface_sharpness =
    this->simulation_parameters.multiphysics.vof_parameters.sharpening
      .interface_sharpness;

  FEValues<dim> fe_values_vof(*this->mapping,
                              *this->fe,
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

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
              // phase_value      = std::min(std::max(phase_value, 0.), 1.);

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

                  // $$ (if 0 <= \phi <= c)  {\Phi = c ^ (1 - \alpha) * (\phi
                  // ^ \alpha)}$$
                  // $$ (if c <  \phi <= 1)  {\Phi = 1 - (1 - c) ^ (1 -
                  // \alpha)
                  // * (1 - \phi) ^ \alpha}
                  if (phase_value >= 0 && phase_value <= sharpening_threshold)
                    local_rhs_phase_fraction(i) +=
                      std::pow(sharpening_threshold,
                               (1. - interface_sharpness)) *
                      std::pow(phase_value, interface_sharpness) *
                      phi_phase[i] * fe_values_vof.JxW(q);
                  else
                    {
                      local_rhs_phase_fraction(i) +=
                        (1 -
                         std::pow((1. - sharpening_threshold),
                                  (1. - interface_sharpness)) *
                           std::pow((1. - phase_value), interface_sharpness)) *
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

  if (this->simulation_parameters.multiphysics.vof_parameters.sharpening
        .verbosity == Parameters::Verbosity::extra_verbose)
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

  if (this->simulation_parameters.multiphysics.vof_parameters.sharpening
        .verbosity == Parameters::Verbosity::extra_verbose)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  this->nonzero_constraints.distribute(
    completely_distributed_phase_fraction_solution);
  solution = completely_distributed_phase_fraction_solution;
}

template <int dim>
void
VolumeOfFluid<dim>::assemble_mass_matrix(
  TrilinosWrappers::SparseMatrix &mass_matrix)
{
  QGauss<dim> quadrature_formula(this->cell_quadrature->size());

  FEValues<dim> fe_values_vof(*this->fe,
                              quadrature_formula,
                              update_values | update_JxW_values);


  const unsigned int dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_vof.reinit(cell);
          cell_matrix = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
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
  this->nb_cells_wet    = 0;
  this->nb_cells_peeled = 0;

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
              if (this->simulation_parameters.multiphysics
                    .use_time_average_velocity_field &&
                  simulation_control->get_current_time() >
                    this->simulation_parameters.post_processing.initial_time)
                {
                  apply_peeling_wetting(
                    i_bc,
                    *multiphysics->get_block_time_average_solution(
                      PhysicsID::fluid_dynamics));
                }
              else
                {
                  apply_peeling_wetting(i_bc,
                                        *multiphysics->get_block_solution(
                                          PhysicsID::fluid_dynamics));
                }
            }
          else
            {
              if (this->simulation_parameters.multiphysics
                    .use_time_average_velocity_field &&
                  simulation_control->get_current_time() >
                    this->simulation_parameters.post_processing.initial_time)
                {
                  apply_peeling_wetting(
                    i_bc,
                    *multiphysics->get_time_average_solution(
                      PhysicsID::fluid_dynamics));
                }
              else
                {
                  apply_peeling_wetting(i_bc,
                                        *multiphysics->get_solution(
                                          PhysicsID::fluid_dynamics));
                }
            }
        }
    } // end loop on boundary_conditions_vof

  // Output total of peeled/wet cells in the entire domain
  if (this->simulation_parameters.multiphysics.vof_parameters.peeling_wetting
        .verbosity != Parameters::Verbosity::quiet)
    {
      auto mpi_communicator = this->triangulation->get_communicator();

      this->pcout << "Peeling/wetting correction at step "
                  << this->simulation_control->get_step_number() << std::endl;
      this->pcout << "  -number of wet cells: "
                  << Utilities::MPI::sum(this->nb_cells_wet, mpi_communicator)
                  << std::endl;
      this->pcout << "  -number of peeled cells: "
                  << Utilities::MPI::sum(this->nb_cells_peeled,
                                         mpi_communicator)
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
                                       update_values);

  FEFaceValues<dim> fe_face_values_fd(*this->mapping,
                                      dof_handler_fd->get_fe(),
                                      *this->face_quadrature,
                                      update_values | update_gradients |
                                        update_normal_vectors);

  const unsigned int n_q_points_face = this->face_quadrature->size();

  const FEValuesExtractors::Scalar pressure(dim);
  std::vector<double>              pressure_values(n_q_points_face);
  std::vector<Tensor<1, dim>>      pressure_gradients(n_q_points_face);
  std::vector<double>              phase_values(n_q_points_face);
  Tensor<1, dim>                   normal_vector_fd;

  unsigned int boundary_id =
    this->simulation_parameters.boundary_conditions_vof.id[i_bc];

  // Physical properties
  const auto density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::map<field, std::vector<double>> fields;

  std::vector<double> density_0(n_q_points_face);
  std::vector<double> density_1(n_q_points_face);

  // Parse fluid present solution to calculate average pressure
  double pressure_monitored_avg = 0.;

  if (multiphysics->fluid_dynamics_is_block())
    {
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          pressure_monitored_avg = find_monitored_fluid_avg_pressure(
            this->present_solution,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.multiphysics.vof_parameters.conservation
              .monitored_fluid);
        }
      else
        {
          pressure_monitored_avg = find_monitored_fluid_avg_pressure(
            this->present_solution,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.multiphysics.vof_parameters.conservation
              .monitored_fluid);
        }
    }
  else
    {
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          pressure_monitored_avg = find_monitored_fluid_avg_pressure(
            this->present_solution,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.multiphysics.vof_parameters.conservation
              .monitored_fluid);
        }
      else
        {
          pressure_monitored_avg = find_monitored_fluid_avg_pressure(
            this->present_solution,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.multiphysics.vof_parameters.conservation
              .monitored_fluid);
        }
    }

  // Loop on cell_vof
  for (const auto &cell_vof : dof_handler.active_cell_iterators())
    {
      if (cell_vof->is_locally_owned() && cell_vof->at_boundary())
        {
          // Initialize sum of values on quadrature points
          unsigned int phase_denser_fluid_q(0);
          unsigned int phase_lighter_fluid_q(0);
          float        phase_values_q(0);
          double       pressure_values_q(0);
          unsigned int n_q_negative_pressure_grad(0);
          unsigned int n_q_positive_pressure_grad(0);

          bool cell_face_at_pw_bounday(false);

          // Local index (see deal.II step 4)
          cell_vof->get_dof_indices(dof_indices_vof);

          for (const auto face : cell_vof->face_indices())
            {
              if (cell_vof->face(face)->at_boundary() &&
                  cell_vof->face(face)->boundary_id() == boundary_id)
                {
                  cell_face_at_pw_bounday = true;
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
                  for (unsigned int q = 0; q < n_q_points_face; q++)
                    {
                      // Get denser/lighter fluid id
                      if (density_1[q] > density_0[q])
                        phase_denser_fluid_q += 1;
                      else
                        phase_lighter_fluid_q += 1;

                      // Condition on the pressure gradient
                      bool peeling_condition_q(false);
                      bool wetting_condition_q(false);
                      for (unsigned int d = 0; d < dim; d++)
                        {
                          double pressure_gradient_q =
                            pressure_gradients[q][d] *
                            fe_face_values_fd.normal_vector(q)[d];

                          // Peeling if pressure_gradient_q is negative
                          if (pressure_gradient_q <
                              -std::numeric_limits<double>::min())
                            peeling_condition_q = true;
                          // Wetting if pressure_gradient_q is positive
                          if (pressure_gradient_q >
                              std::numeric_limits<double>::min())
                            wetting_condition_q = true;
                        }

                      // Check if condition is reached for this quadrature point
                      if (peeling_condition_q)
                        n_q_negative_pressure_grad += 1;
                      if (wetting_condition_q)
                        n_q_positive_pressure_grad += 1;

                      pressure_values_q += pressure_values[q];
                      phase_values_q += phase_values[q];
                    } // end loop on quadrature points
                }     // end condition face at pw boundary
            }         // end loop on faces

          // Check if cell is to be treated
          if (cell_face_at_pw_bounday)
            {
              // Caculate average values on the cell
              unsigned int phase_denser_fluid_cell =
                round(phase_denser_fluid_q / n_q_points_face);
              unsigned int phase_lighter_fluid_cell =
                round(phase_lighter_fluid_q / n_q_points_face);
              double phase_values_cell =
                phase_values_q / static_cast<double>(n_q_points_face);
              double pressure_values_cell =
                pressure_values_q / static_cast<double>(n_q_points_face);

              // Check peeling condition
              // Note that it has priority over the wetting condition
              if (pressure_values_cell < pressure_monitored_avg and
                  n_q_negative_pressure_grad > n_q_points_face * 0.5)
                {
                  if (this->simulation_parameters.multiphysics.vof_parameters
                        .peeling_wetting.enable_peeling)
                    {
                      // Peel the higher density fluid
                      if (phase_denser_fluid_cell == 1 and
                          phase_values_cell > 0.5)
                        {
                          // Progressive phase value change
                          double new_phase =
                            std::max(phase_values_cell - 0.25,
                                     static_cast<double>(
                                       phase_lighter_fluid_cell));

                          change_cell_phase(PhaseChange::peeling,
                                            new_phase,
                                            solution_pw,
                                            dof_indices_vof);
                        }
                      if (phase_denser_fluid_cell == 0 and
                          phase_values_cell < 0.5)
                        {
                          // Progressive phase value change
                          double new_phase =
                            std::max(phase_values_cell + 0.25,
                                     static_cast<double>(
                                       phase_lighter_fluid_cell));

                          change_cell_phase(PhaseChange::peeling,
                                            new_phase,
                                            solution_pw,
                                            dof_indices_vof);
                        }
                    }
                } // end condition peeling
              // Check wetting condition if the cell is not marked as pelt
              else if (pressure_values_cell > pressure_monitored_avg and
                       n_q_positive_pressure_grad > n_q_points_face * 0.5)
                {
                  if (this->simulation_parameters.multiphysics.vof_parameters
                        .peeling_wetting.enable_wetting)
                    {
                      if (phase_denser_fluid_cell == 1 and
                          phase_values_cell > 0.5)
                        {
                          // Progressive phase value change
                          double new_phase =
                            std::min(phase_values_cell + 0.25,
                                     static_cast<double>(
                                       phase_denser_fluid_cell));

                          change_cell_phase(PhaseChange::wetting,
                                            new_phase,
                                            solution_pw,
                                            dof_indices_vof);
                        }
                      else if (phase_denser_fluid_cell == 0 and
                               phase_values_cell < 0.5)
                        {
                          // Progressive phase value change
                          double new_phase =
                            std::min(phase_values_cell - 0.25,
                                     static_cast<double>(
                                       phase_denser_fluid_cell));

                          change_cell_phase(PhaseChange::wetting,
                                            new_phase,
                                            solution_pw,
                                            dof_indices_vof);
                        }
                    }
                } // end condition wetting condition
            }     // end condition cell_face_at_pw_bounday

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
  const double &                              new_phase,
  TrilinosWrappers::MPI::Vector &             solution_pw,
  const std::vector<types::global_dof_index> &dof_indices_vof)
{
  if (type == PhaseChange::wetting)
    {
      for (unsigned int k = 0; k < fe->dofs_per_cell; ++k)
        {
          solution_pw[dof_indices_vof[k]]     = new_phase;
          this->marker_pw[dof_indices_vof[k]] = 1;
        }
      this->nb_cells_wet++;
    }
  else if (type == PhaseChange::peeling)
    {
      for (unsigned int k = 0; k < fe->dofs_per_cell; ++k)
        {
          solution_pw[dof_indices_vof[k]]     = new_phase;
          this->marker_pw[dof_indices_vof[k]] = -1;
        }
      this->nb_cells_peeled++;
    }
}

template <int dim>
void
VolumeOfFluid<dim>::apply_phase_filter()
{
  // Initializations
  auto mpi_communicator = this->triangulation->get_communicator();
  TrilinosWrappers::MPI::Vector filtered_solution_owned(
    this->locally_owned_dofs, mpi_communicator);
  filtered_solution_owned = this->present_solution;
  filtered_solution.reinit(this->present_solution);

  // Create filter object
  filter = VolumeOfFluidFilterBase::model_cast(
    this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

  // Apply filter to solution
  for (auto p : this->locally_owned_dofs)
    {
      filtered_solution_owned[p] =
        filter->filter_phase(filtered_solution_owned[p]);
    }
  filtered_solution = filtered_solution_owned;

  if (this->simulation_parameters.multiphysics.vof_parameters.phase_filter
        .verbosity == Parameters::Verbosity::verbose)
    {
      this->pcout << "Filtered phase values: " << std::endl;
      for (const double filtered_phase : filtered_solution)
        {
          this->pcout << filtered_phase << std::endl;
        }
    }
}

template class VolumeOfFluid<2>;
template class VolumeOfFluid<3>;
