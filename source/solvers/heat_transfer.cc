#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/heat_transfer.h>
#include <solvers/heat_transfer_assemblers.h>
#include <solvers/heat_transfer_scratch_data.h>
#include <solvers/postprocessing_cfd.h>

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

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
void
HeatTransfer<dim>::assemble_matrix_and_rhs()
{
  assemble_system_matrix();
  assemble_system_rhs();

  if (this->simulation_parameters.nitsche->number_solids > 0)
    assemble_nitsche_heat_restriction(true);
}


template <int dim>
void
HeatTransfer<dim>::assemble_rhs()
{
  assemble_system_rhs();
  if (this->simulation_parameters.nitsche->number_solids > 0)
    assemble_nitsche_heat_restriction(true);
}

template <int dim>
void
HeatTransfer<dim>::assemble_nitsche_heat_restriction(bool assemble_matrix)
{
  Assert(
    !this->simulation_parameters.physical_properties_manager.is_non_newtonian(),
    RequiresConstantViscosity("assemble_nitsche_heat_restriction"));

  // Evaluate fluid properties
  // auto density_model =
  //  this->simulation_parameters.physical_properties_manager.get_density();
  // auto specific_heat_model =
  //  this->simulation_parameters.physical_properties_manager.get_specific_heat();
  // auto conductivity_model =
  //  this->simulation_parameters.physical_properties_manager
  //    .get_thermal_conductivity();
  // std::map<field, double> field_values;

  // double rho_cp = density_model->value(field_values) *
  //                 specific_heat_model->value(field_values);

  auto solids = *this->multiphysics->get_solids(
    this->simulation_parameters.nitsche->number_solids);

  // Loops over solids
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      // Initialize the solid and its particles
      std::shared_ptr<Particles::ParticleHandler<dim>> &solid_ph =
        solids[i_solid]->get_solid_particle_handler();

      if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
            ->enable_heat_bc)
        {
          const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

          std::vector<types::global_dof_index> heat_dof_indices(dofs_per_cell);
          FullMatrix<double>     local_matrix(dofs_per_cell, dofs_per_cell);
          dealii::Vector<double> local_rhs(dofs_per_cell);

          Function<dim> *solid_temperature =
            solids[i_solid]->get_solid_temperature();

          // Penalization terms
          const double beta =
            this->simulation_parameters.nitsche->nitsche_solids[i_solid]
              ->beta_heat;

          // Loop over all local particles
          auto particle = solid_ph->begin();
          while (particle != solid_ph->end())
            {
              local_matrix = 0;
              local_rhs    = 0;

              const auto &cell = particle->get_surrounding_cell();

              // double h_cell = 0;
              // if (dim == 2)
              //   h_cell =
              //     std::sqrt(4. * cell->measure() / M_PI) /
              //     this->simulation_parameters.fem_parameters.velocity_order;
              // else if (dim == 3)
              //   h_cell =
              //     pow(6 * cell->measure() / M_PI, 1. / 3.) /
              //     this->simulation_parameters.fem_parameters.velocity_order;
              //  Penalty parameter is disabled for heat transfer from since we
              //  don't necessarily want to strictly impose a Dirichlet BC
              const double penalty_parameter =
                1.; /// std::pow(h_cell * h_cell, double(dim) / double(dim));
              const auto &dh_cell =
                typename DoFHandler<dim>::cell_iterator(*cell,
                                                        &this->dof_handler);
              dh_cell->get_dof_indices(heat_dof_indices);

              const auto pic = solid_ph->particles_in_cell(cell);
              Assert(pic.begin() == particle, ExcInternalError());
              for (const auto &p : pic)
                {
                  const ArrayView<const double> properties = p.get_properties();

                  double      fluid_temperature = 0;
                  const auto &ref_q             = p.get_reference_location();
                  const auto &real_q            = p.get_location();
                  const auto &JxW               = properties[0];

                  for (unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                      // Get the fluid temperature at non-quadrature point
                      // (particle in fluid)
                      auto &evaluation_point = this->evaluation_point;
                      fluid_temperature +=
                        evaluation_point[heat_dof_indices[k]] *
                        this->fe->shape_value(k, ref_q);
                    }
                  // Assemble the matrix or the rhs
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      if (assemble_matrix)
                        {
                          for (unsigned int j = 0; j < dofs_per_cell; ++j)
                            {
                              local_matrix(i, j) +=
                                penalty_parameter * beta *
                                this->fe->shape_value(i, ref_q) *
                                this->fe->shape_value(j, ref_q) * JxW;
                            }
                        }
                      else
                        {
                          // Regular residual
                          local_rhs(i) += penalty_parameter * beta *
                                          (solid_temperature->value(real_q, 0) -
                                           fluid_temperature) *
                                          this->fe->shape_value(i, ref_q) * JxW;
                        }
                    }
                }
              const AffineConstraints<double> &constraints_used =
                this->zero_constraints;
              auto &system_rhs = this->system_rhs;
              constraints_used.distribute_local_to_global(local_matrix,
                                                          local_rhs,
                                                          heat_dof_indices,
                                                          this->system_matrix,
                                                          system_rhs);
              particle = pic.end();
            }
          this->system_matrix.compress(VectorOperation::add);
          this->system_rhs.compress(VectorOperation::add);
        }
    }
}

template <int dim>
void
HeatTransfer<dim>::postprocess_heat_flux_on_nitsche_ib()
{
  auto solids = *this->multiphysics->get_solids(
    this->simulation_parameters.nitsche->number_solids);
  std::vector<double> heat_flux_on_nitsche_ib_vector(solids.size(), 0);

  // Loops over solids
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      // Initialize the solid and its particles
      std::shared_ptr<Particles::ParticleHandler<dim>> &solid_ph =
        solids[i_solid]->get_solid_particle_handler();
      double heat_flux_on_nitsche_bc = 0;

      if (this->simulation_parameters.nitsche->nitsche_solids[i_solid]
            ->enable_heat_bc)
        {
          const unsigned int dofs_per_cell = this->fe->dofs_per_cell;

          std::vector<types::global_dof_index> heat_dof_indices(dofs_per_cell);
          FullMatrix<double>     local_matrix(dofs_per_cell, dofs_per_cell);
          dealii::Vector<double> local_rhs(dofs_per_cell);

          Function<dim> *solid_temperature =
            solids[i_solid]->get_solid_temperature();

          // Penalization terms
          const double beta =
            this->simulation_parameters.nitsche->nitsche_solids[i_solid]
              ->beta_heat;

          // Loop over all local particles
          auto particle = solid_ph->begin();
          while (particle != solid_ph->end())
            {
              local_matrix = 0;
              local_rhs    = 0;

              const auto &cell = particle->get_surrounding_cell();

              const auto &dh_cell =
                typename DoFHandler<dim>::cell_iterator(*cell,
                                                        &this->dof_handler);
              dh_cell->get_dof_indices(heat_dof_indices);

              const auto pic = solid_ph->particles_in_cell(cell);
              Assert(pic.begin() == particle, ExcInternalError());
              for (const auto &p : pic)
                {
                  const ArrayView<const double> properties = p.get_properties();

                  double      fluid_temperature = 0;
                  const auto &ref_q             = p.get_reference_location();
                  const auto &real_q            = p.get_location();
                  const auto &JxW               = properties[0];

                  for (unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                      // Get the fluid temperature at non-quadrature point
                      // (particle in fluid)
                      fluid_temperature +=
                        this->present_solution[heat_dof_indices[k]] *
                        this->fe->shape_value(k, ref_q);
                    }
                  // Calculate heat flux on the Nitsche IB
                  heat_flux_on_nitsche_bc +=
                    beta *
                    (solid_temperature->value(real_q, 0) - fluid_temperature) *
                    JxW;
                }
              particle = pic.end();
            }
        }
      heat_flux_on_nitsche_ib_vector[i_solid] =
        Utilities::MPI::sum(heat_flux_on_nitsche_bc,
                            triangulation->get_communicator());
    }

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Heat flux at the nitsche immersed boundary : "
                  << std::endl;
      for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
        {
          this->pcout << "\t Nitsche boundary " << i_solid << " : "
                      << heat_flux_on_nitsche_ib_vector[i_solid] << std::endl;
        }
    }

  // Filling table
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      this->heat_flux_table.add_value("flux_nitsche_solid_" +
                                        Utilities::int_to_string(i_solid, 1),
                                      heat_flux_on_nitsche_ib_vector[i_solid]);
    }
}

template <int dim>
void
HeatTransfer<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Laser heat source
  if (this->simulation_parameters.laser_parameters->activate_laser)
    {
      if (this->simulation_parameters.multiphysics.VOF)
        {
          // Call for the specific assembler
          // Assembler of the laser source term applied only to the metal phase
          this->assemblers.push_back(
            std::make_shared<HeatTransferAssemblerLaserVOF<dim>>(
              this->simulation_control,
              this->simulation_parameters.laser_parameters));

          // Assembler of the radiation sink term applied only at the air/metal
          // interface. The radiation term in that case is treated as a source
          // term instead of a boundary term.
          if (this->simulation_parameters.laser_parameters->radiation
                .enable_radiation)
            {
              this->assemblers.push_back(
                std::make_shared<
                  HeatTransferAssemblerFreeSurfaceRadiationVOF<dim>>(
                  this->simulation_control,
                  this->simulation_parameters.laser_parameters));
            }
        }
      else
        {
          this->assemblers.push_back(
            std::make_shared<HeatTransferAssemblerLaser<dim>>(
              this->simulation_control,
              this->simulation_parameters.laser_parameters));
        }
    }

  // Robin boundary condition
  this->assemblers.push_back(
    std::make_shared<HeatTransferAssemblerRobinBC<dim>>(
      this->simulation_control, simulation_parameters.boundary_conditions_ht));

  if (this->simulation_parameters.multiphysics.viscous_dissipation)
    {
      if (this->simulation_parameters.multiphysics.VOF)
        {
          // Call for the specific assembler
          this->assemblers.push_back(
            std::make_shared<HeatTransferAssemblerViscousDissipationVOF<dim>>(
              this->simulation_control,
              this->simulation_parameters.multiphysics.vof_parameters
                .viscous_dissipative_fluid));
        }
      else
        {
          this->assemblers.push_back(
            std::make_shared<HeatTransferAssemblerViscousDissipation<dim>>(
              this->simulation_control));
        }
    }

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(
        std::make_shared<HeatTransferAssemblerBDF<dim>>(
          this->simulation_control));
    }

  // Core assembler
  this->assemblers.push_back(
    std::make_shared<HeatTransferAssemblerCore<dim>>(this->simulation_control));
}

template <int dim>
void
HeatTransfer<dim>::assemble_system_matrix()
{
  this->system_matrix = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = HeatTransferScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->temperature_mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->temperature_mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);
    }

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &HeatTransfer::assemble_local_system_matrix,
                  &HeatTransfer::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);

  if (this->simulation_parameters.nitsche->number_solids > 0)
    assemble_nitsche_heat_restriction(true);
}

template <int dim>
void
HeatTransfer<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  HeatTransferScratchData<dim> &                        scratch_data,
  StabilizedMethodsCopyData &                           copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  auto &source_term = simulation_parameters.source_term->heat_transfer_source;
  source_term.set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      &source_term);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

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
              PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(velocity_cell,
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
          scratch_data.reinit_velocity(velocity_cell,
                                       *multiphysics->get_time_average_solution(
                                         PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(velocity_cell,
                                       *multiphysics->get_solution(
                                         PhysicsID::fluid_dynamics));
        }
    }

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(phase_cell,
                              *this->multiphysics->get_filtered_solution(
                                PhysicsID::VOF),
                              std::vector<TrilinosWrappers::MPI::Vector>());
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
HeatTransfer<dim>::copy_local_matrix_to_global_matrix(
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
HeatTransfer<dim>::assemble_system_rhs()
{
  // TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = HeatTransferScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->temperature_mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->temperature_mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);
    }

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &HeatTransfer::assemble_local_system_rhs,
                  &HeatTransfer::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_parameters.nitsche->number_solids > 0)
    assemble_nitsche_heat_restriction(false);
}

template <int dim>
void
HeatTransfer<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  HeatTransferScratchData<dim> &                        scratch_data,
  StabilizedMethodsCopyData &                           copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  auto &source_term = simulation_parameters.source_term->heat_transfer_source;
  source_term.set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      this->solution_stages,
                      &source_term);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

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
              PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(velocity_cell,
                                       *multiphysics->get_block_solution(
                                         PhysicsID::fluid_dynamics));
        }
      scratch_data.reinit_velocity_gradient(
        *multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
    }
  else
    {
      if (this->simulation_parameters.multiphysics
            .use_time_average_velocity_field &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing.initial_time)
        {
          scratch_data.reinit_velocity(velocity_cell,
                                       *multiphysics->get_time_average_solution(
                                         PhysicsID::fluid_dynamics));
        }
      else
        {
          scratch_data.reinit_velocity(velocity_cell,
                                       *multiphysics->get_solution(
                                         PhysicsID::fluid_dynamics));
        }

      scratch_data.reinit_velocity_gradient(
        *multiphysics->get_solution(PhysicsID::fluid_dynamics));
    }

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(phase_cell,
                              *this->multiphysics->get_filtered_solution(
                                PhysicsID::VOF),
                              std::vector<TrilinosWrappers::MPI::Vector>());
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
HeatTransfer<dim>::copy_local_rhs_to_global_rhs(
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
HeatTransfer<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  data_out.add_data_vector(dof_handler, present_solution, "temperature");
}

template <int dim>
double
HeatTransfer<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values(*this->temperature_mapping,
                          *fe,
                          *this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = this->cell_quadrature->size();

  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->temperature;
  exact_solution.set_time(simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_scalar_values);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values.get_quadrature_points(),
                                    q_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = q_scalar_values[q];
              double exact = q_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
void
HeatTransfer<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity !=
        Parameters::Verbosity::quiet)
    {
      error_table.omit_column_from_convergence_rate_evaluation("cells");


      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_temperature", true);
      error_table.set_precision("error_temperature",
                                simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
HeatTransfer<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
HeatTransfer<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      double temperature_error = calculate_L2_error();

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());
      error_table.add_value("error_temperature", temperature_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error temperature : " << temperature_error
                      << std::endl;
        }
    }

  // Set-up domain name for output files
  Parameters::FluidIndicator monitored_fluid =
    this->simulation_parameters.post_processing.postprocessed_fluid;
  // default: monophase simulations
  std::string domain_name("fluid");
  bool        gather_vof(false);

  if (this->simulation_parameters.physical_properties_manager
        .get_number_of_fluids() == 2)
    {
      // Multiphase flow
      gather_vof = true;
      switch (monitored_fluid)
        {
            default: {
              domain_name = "fluid";
              break;
            }
            case Parameters::FluidIndicator::fluid0: {
              domain_name = "fluid_0";
              break;
            }
            case Parameters::FluidIndicator::fluid1: {
              domain_name = "fluid_1";
              break;
            }
        }
    }

  // Temperature statistics
  if (simulation_parameters.post_processing.calculate_temperature_statistics)
    {
      postprocess_temperature_statistics(gather_vof,
                                         monitored_fluid,
                                         domain_name);

      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_temperature_statistics(domain_name);
    }

  // Heat flux calculation
  if (simulation_parameters.post_processing.calculate_heat_flux)
    {
      // Parse fluid present solution
      if (multiphysics->fluid_dynamics_is_block())
        {
          postprocess_heat_flux_on_bc(gather_vof,
                                      *multiphysics->get_block_solution(
                                        PhysicsID::fluid_dynamics));

          postprocess_thermal_energy_in_fluid(gather_vof,
                                              monitored_fluid,
                                              domain_name,
                                              *multiphysics->get_block_solution(
                                                PhysicsID::fluid_dynamics));
        }
      else
        {
          postprocess_heat_flux_on_bc(
            gather_vof, *multiphysics->get_solution(PhysicsID::fluid_dynamics));

          postprocess_thermal_energy_in_fluid(gather_vof,
                                              monitored_fluid,
                                              domain_name,
                                              *multiphysics->get_solution(
                                                PhysicsID::fluid_dynamics));
        }

      if (this->simulation_parameters.nitsche->number_solids > 0)
        postprocess_heat_flux_on_nitsche_ib();

      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_heat_flux(domain_name);
    }
}

template <int dim>
void
HeatTransfer<dim>::pre_mesh_adaptation()
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
HeatTransfer<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_communicator();


  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer->interpolate(tmp);

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
HeatTransfer<dim>::compute_kelly(
  const std::pair<const Parameters::MeshAdaptation::Variable,
                  Parameters::MultipleAdaptationParameters> &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Parameters::MeshAdaptation::Variable::temperature)
    {
      const FEValuesExtractors::Scalar temperature(0);

      KellyErrorEstimator<dim>::estimate(
        *this->temperature_mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(temperature));
    }
}

template <int dim>
void
HeatTransfer<dim>::write_checkpoint()
{
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;

  solution_transfer =
    std::make_shared<parallel::distributed::
                       SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>>(
      dof_handler);

  sol_set_transfer.push_back(&present_solution);
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&previous_solutions[i]);
    }
  solution_transfer->prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
HeatTransfer<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading heat transfer checkpoint" << std::endl;

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

  solution_transfer->deserialize(input_vectors);

  present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }
}


template <int dim>
void
HeatTransfer<dim>::setup_dofs()
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

    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions_ht.size;
         ++i_bc)
      {
        // Dirichlet condition : imposed temperature at i_bc
        if (this->simulation_parameters.boundary_conditions_ht.type[i_bc] ==
            BoundaryConditions::BoundaryType::temperature)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_ht.id[i_bc],
              dealii::Functions::ConstantFunction<dim>(
                this->simulation_parameters.boundary_conditions_ht.value[i_bc]),
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

    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions_ht.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions_ht.type[i_bc] ==
            BoundaryConditions::BoundaryType::temperature)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_ht.id[i_bc],
              Functions::ZeroFunction<dim>(),
              zero_constraints);
          }
      }
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
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of thermal degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the heat transfer dof_handler and present solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::heat_transfer, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::heat_transfer, &this->present_solution);
  multiphysics->set_previous_solutions(PhysicsID::heat_transfer,
                                       &this->previous_solutions);
}

template <int dim>
void
HeatTransfer<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*this->temperature_mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->temperature,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  percolate_time_vectors();
}

template <int dim>
void
HeatTransfer<dim>::solve_linear_system(const bool initial_step,
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

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template <int dim>
void
HeatTransfer<dim>::postprocess_temperature_statistics(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double> local_temperature_values(n_q_points);
  FEValues<dim>       fe_values_ht(*this->temperature_mapping,
                             *this->fe,
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Initialize VOF information
  const DoFHandler<dim> *        dof_handler_vof = NULL;
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  std::vector<double>            filtered_phase_values(n_q_points);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->temperature_mapping,
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }

  // Other initializations
  double phase_coefficient(0.);
  double volume_integral(0.);
  double temperature_integral(0.);
  double minimum_temperature(std::numeric_limits<double>::max());
  double maximum_temperature(std::numeric_limits<double>::lowest());
  // boolean used to state if the quadrature point is in the
  // postprocessed_fluid, i.e. its temperature should be considered
  // in the minimum and maximum temperature calculation
  bool point_is_in_postprocessed_fluid(false);

  // Calculate min, max and average
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Gather heat transfer information
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              std::tie(phase_coefficient, point_is_in_postprocessed_fluid) =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q]);

              if (point_is_in_postprocessed_fluid)
                {
                  maximum_temperature =
                    std::max(local_temperature_values[q], maximum_temperature);
                  minimum_temperature =
                    std::min(local_temperature_values[q], minimum_temperature);
                }

              volume_integral += fe_values_ht.JxW(q) * phase_coefficient;
              temperature_integral += local_temperature_values[q] *
                                      phase_coefficient * fe_values_ht.JxW(q);
            } // end loop on quadrature points
        }
    } // end loop on cell

  minimum_temperature =
    Utilities::MPI::min(minimum_temperature, mpi_communicator);
  maximum_temperature =
    Utilities::MPI::max(maximum_temperature, mpi_communicator);

  volume_integral = Utilities::MPI::sum(volume_integral, mpi_communicator);
  temperature_integral =
    Utilities::MPI::sum(temperature_integral, mpi_communicator);
  double temperature_average = temperature_integral / volume_integral;

  // Calculate standard deviation
  double temperature_variance_integral = 0;
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              phase_coefficient =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q])
                  .first;

              temperature_variance_integral +=
                (local_temperature_values[q] - temperature_average) *
                (local_temperature_values[q] - temperature_average) *
                phase_coefficient * fe_values_ht.JxW(q);
            } // end loop on quadrature points
        }
    } // end loop on cell

  temperature_variance_integral =
    Utilities::MPI::sum(temperature_variance_integral, mpi_communicator);
  double temperature_variance = temperature_variance_integral / volume_integral;
  double temperature_std_deviation = std::sqrt(temperature_variance);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Temperature statistics on " << domain_name << " : "
                  << std::endl;
      this->pcout << "\t     Min : " << minimum_temperature << std::endl;
      this->pcout << "\t     Max : " << maximum_temperature << std::endl;
      this->pcout << "\t Average : " << temperature_average << std::endl;
      this->pcout << "\t Std-Dev : " << temperature_std_deviation << std::endl;
    }

  // Filling table
  this->statistics_table.add_value(
    "time", this->simulation_control->get_current_time());
  this->statistics_table.add_value("min", minimum_temperature);
  this->statistics_table.add_value("max", maximum_temperature);
  this->statistics_table.add_value("average", temperature_average);
  this->statistics_table.add_value("std-dev", temperature_std_deviation);
}

template <int dim>
void
HeatTransfer<dim>::write_temperature_statistics(const std::string domain_name)
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.temperature_output_name + "_" +
        domain_name + ".dat";
      std::ofstream output(filename.c_str());

      this->statistics_table.write_text(output);
    }
}

template <int dim>
template <typename VectorType>
void
HeatTransfer<dim>::postprocess_heat_flux_on_bc(
  const bool        gather_vof,
  const VectorType &current_solution_fd)
{
  const unsigned int n_q_points_face  = this->face_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double>         local_temperature_values(n_q_points_face);
  std::vector<Tensor<1, dim>> temperature_gradient(n_q_points_face);
  FEFaceValues<dim>           fe_face_values_ht(*this->temperature_mapping,
                                      this->dof_handler.get_fe(),
                                      *this->face_quadrature,
                                      update_values | update_gradients |
                                        update_JxW_values |
                                        update_normal_vectors);

  // Initialize fluid dynamics information
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      local_velocity_values(n_q_points_face);
  FEFaceValues<dim>                fe_face_values_fd(*this->temperature_mapping,
                                      dof_handler_fd->get_fe(),
                                      *this->face_quadrature,
                                      update_values);

  // Initialize VOF information
  DoFHandler<dim> *                  dof_handler_vof = NULL;
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_vof;
  std::vector<double>                filtered_phase_values(n_q_points_face);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_face_values_vof =
        std::make_shared<FEFaceValues<dim>>(*this->temperature_mapping,
                                            dof_handler_vof->get_fe(),
                                            *this->face_quadrature,
                                            update_values);
    }

  // Initialize fluid properties
  auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::map<field, double>              field_values;
  std::map<field, std::vector<double>> fields;

  // monophase flow
  double density(0.);
  double specific_heat(0.);
  double thermal_conductivity(0.);

  // multiphase flow
  std::vector<double> density_0(n_q_points_face);
  std::vector<double> specific_heat_0(n_q_points_face);
  std::vector<double> thermal_conductivity_0(n_q_points_face);
  std::vector<double> density_1(n_q_points_face);
  std::vector<double> specific_heat_1(n_q_points_face);
  std::vector<double> thermal_conductivity_1(n_q_points_face);

  switch (properties_manager.get_number_of_fluids())
    {
        default: {
          // Get values for monophase flow
          const auto density_model = properties_manager.get_density();
          const auto specific_heat_model =
            properties_manager.get_specific_heat();
          const auto conductivity_model =
            properties_manager.get_thermal_conductivity();

          density              = density_model->value(field_values);
          specific_heat        = specific_heat_model->value(field_values);
          thermal_conductivity = conductivity_model->value(field_values);

          break;
        }
        case 2: {
          // Get prm values for multiphase flow - will be blended in the
          // integration loop
          const auto density_models = properties_manager.get_density_vector();
          const auto specific_heat_models =
            properties_manager.get_specific_heat_vector();
          const auto thermal_conductivity_models =
            properties_manager.get_thermal_conductivity_vector();

          density_models[0]->vector_value(fields, density_0);
          specific_heat_models[0]->vector_value(fields, specific_heat_0);
          thermal_conductivity_models[0]->vector_value(fields,
                                                       thermal_conductivity_0);

          density_models[1]->vector_value(fields, density_1);
          specific_heat_models[1]->vector_value(fields, specific_heat_1);
          thermal_conductivity_models[1]->vector_value(fields,
                                                       thermal_conductivity_1);

          break;
        }
    } // end switch on number of fluids

  std::vector<double> heat_flux_vector(
    this->simulation_parameters.boundary_conditions_ht.size, 0);

  std::vector<double> convective_flux_vector(
    this->simulation_parameters.boundary_conditions_ht.size, 0);


  // Get vector of heat transfer boundary conditions
  const auto boundary_conditions_ht_ids =
    this->simulation_parameters.boundary_conditions_ht.id;

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (const auto face : cell->face_indices())
            {
              if (cell->face(face)->at_boundary())
                {
                  const auto boundary_id =
                    std::find(begin(boundary_conditions_ht_ids),
                              end(boundary_conditions_ht_ids),
                              cell->face(face)->boundary_id());

                  if (boundary_id != end(boundary_conditions_ht_ids))
                    {
                      unsigned int vector_index =
                        boundary_id - boundary_conditions_ht_ids.begin();

                      // Gather h coefficient and T_inf
                      const double h_coefficient =
                        this->simulation_parameters.boundary_conditions_ht
                          .h[vector_index];
                      const double T_inf =
                        this->simulation_parameters.boundary_conditions_ht
                          .Tinf[vector_index];


                      // Gather heat transfer information
                      fe_face_values_ht.reinit(cell, face);
                      fe_face_values_ht.get_function_values(
                        this->present_solution, local_temperature_values);
                      fe_face_values_ht.get_function_gradients(
                        this->present_solution, temperature_gradient);

                      // Get fluid dynamics active cell iterator
                      typename DoFHandler<dim>::active_cell_iterator cell_fd(
                        &(*(this->triangulation)),
                        cell->level(),
                        cell->index(),
                        dof_handler_fd);

                      // Gather fluid dynamics information
                      fe_face_values_fd.reinit(cell_fd, face);
                      fe_face_values_fd[velocities].get_function_values(
                        current_solution_fd, local_velocity_values);

                      if (gather_vof)
                        {
                          // Get VOF active cell iterator
                          typename DoFHandler<dim>::active_cell_iterator
                            cell_vof(&(*(this->triangulation)),
                                     cell->level(),
                                     cell->index(),
                                     dof_handler_vof);

                          // Gather VOF information
                          fe_face_values_vof->reinit(cell_vof, face);
                          fe_face_values_vof->get_function_values(
                            *this->multiphysics->get_filtered_solution(
                              PhysicsID::VOF),
                            filtered_phase_values);
                        }

                      // Loop on the quadrature points
                      for (unsigned int q = 0; q < n_q_points_face; q++)
                        {
                          if (properties_manager.get_number_of_fluids() == 2)
                            {
                              // Blend the physical properties using the VOF
                              // field
                              thermal_conductivity = calculate_point_property(
                                filtered_phase_values[q],
                                thermal_conductivity_0[q],
                                thermal_conductivity_1[q]);

                              density = calculate_point_property(
                                filtered_phase_values[q],
                                density_0[q],
                                density_1[q]);

                              specific_heat = calculate_point_property(
                                filtered_phase_values[q],
                                specific_heat_0[q],
                                specific_heat_1[q]);
                            }

                          Tensor<1, dim> normal_vector_ht =
                            -fe_face_values_ht.normal_vector(q);

                          heat_flux_vector[vector_index] +=
                            (-thermal_conductivity * temperature_gradient[q] *
                               normal_vector_ht +
                             local_temperature_values[q] * density *
                               specific_heat * local_velocity_values[q] *
                               normal_vector_ht) *
                            fe_face_values_ht.JxW(q);

                          convective_flux_vector[vector_index] +=
                            h_coefficient *
                            (local_temperature_values[q] - T_inf) *
                            fe_face_values_ht.JxW(q);

                        } // end loop on quadrature points
                    }     // end condition face at heat transfer boundary
                }         // end loop on faces
            }             // End face is a boundary face
        }                 // end condition cell at boundary
    }                     // end loop on cells


  // Sum accross all cores
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_ht.size;
       ++i_bc)
    {
      heat_flux_vector[i_bc] =
        Utilities::MPI::sum(heat_flux_vector[i_bc], mpi_communicator);

      convective_flux_vector[i_bc] =
        Utilities::MPI::sum(convective_flux_vector[i_bc], mpi_communicator);
    }

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout
        << "Total heat flux at the heat transfer boundary conditions : "
        << std::endl;
      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions_ht.size;
           ++i_bc)
        this->pcout << "\t boundary " << i_bc << " : " << heat_flux_vector[i_bc]
                    << std::endl;


      this->pcout
        << "Convective heat flux at the heat transfer boundary conditions : "
        << std::endl;
      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions_ht.size;
           ++i_bc)
        this->pcout << "\t boundary " << i_bc << " : "
                    << convective_flux_vector[i_bc] << std::endl;
    }

  // Filling table
  this->heat_flux_table.add_value("time",
                                  this->simulation_control->get_current_time());
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions_ht.size;
       ++i_bc)
    {
      this->heat_flux_table.add_value("total_flux_bc_" +
                                        Utilities::int_to_string(i_bc, 1),
                                      heat_flux_vector[i_bc]);
      this->heat_flux_table.add_value("convective_flux_bc_" +
                                        Utilities::int_to_string(i_bc, 1),
                                      convective_flux_vector[i_bc]);
    }
}

template void
HeatTransfer<2>::postprocess_heat_flux_on_bc<TrilinosWrappers::MPI::Vector>(
  const bool                           gather_vof,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
HeatTransfer<3>::postprocess_heat_flux_on_bc<TrilinosWrappers::MPI::Vector>(
  const bool                           gather_vof,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
HeatTransfer<2>::postprocess_heat_flux_on_bc<
  TrilinosWrappers::MPI::BlockVector>(
  const bool                                gather_vof,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template void
HeatTransfer<3>::postprocess_heat_flux_on_bc<
  TrilinosWrappers::MPI::BlockVector>(
  const bool                                gather_vof,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template <int dim>
template <typename VectorType>
void
HeatTransfer<dim>::postprocess_thermal_energy_in_fluid(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const VectorType &               current_solution_fd)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double>         local_temperature_values(n_q_points);
  std::vector<Tensor<1, dim>> temperature_gradient(n_q_points);
  FEValues<dim>               fe_values_ht(*this->temperature_mapping,
                             this->dof_handler.get_fe(),
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Initialize fluid dynamics information
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      local_velocity_values(n_q_points);
  FEValues<dim>                    fe_values_fd(*this->temperature_mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  // Initialize VOF information
  const DoFHandler<dim> *        dof_handler_vof = NULL;
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  std::vector<double>            filtered_phase_values(n_q_points);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->temperature_mapping,
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }

  // Initialize fluid properties
  auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::map<field, double>              field_values;
  std::map<field, std::vector<double>> fields;

  // monophase flow
  double density(0.);
  double specific_heat(0.);

  // multiphase flow
  std::vector<double> density_0(n_q_points);
  std::vector<double> specific_heat_0(n_q_points);
  std::vector<double> density_1(n_q_points);
  std::vector<double> specific_heat_1(n_q_points);

  switch (properties_manager.get_number_of_fluids())
    {
        default: {
          // Get values for monophase flow
          const auto density_model = properties_manager.get_density();
          const auto specific_heat_model =
            properties_manager.get_specific_heat();

          density       = density_model->value(field_values);
          specific_heat = specific_heat_model->value(field_values);

          break;
        }
        case 2: {
          // Get prm values for multiphase flow - will be blended in the
          // integration loop
          const auto density_models = properties_manager.get_density_vector();
          const auto specific_heat_models =
            properties_manager.get_specific_heat_vector();

          density_models[0]->vector_value(fields, density_0);
          specific_heat_models[0]->vector_value(fields, specific_heat_0);

          density_models[1]->vector_value(fields, density_1);
          specific_heat_models[1]->vector_value(fields, specific_heat_1);

          break;
        }
    } // end switch on number of fluids

  // Other initializations
  double phase_coefficient(0.);
  double heat_in_domain(0.);

  // Integrate on all domain
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Gather heat transfer information
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_fd);

          // Gather fluid dynamics information
          fe_values_fd.reinit(cell_fd);
          fe_values_fd[velocities].get_function_values(current_solution_fd,
                                                       local_velocity_values);

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          // Loop on the quadrature points
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              phase_coefficient =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q])
                  .first;

              if (properties_manager.get_number_of_fluids() == 2)
                {
                  // Blend the physical properties using the VOF
                  // field
                  density = calculate_point_property(filtered_phase_values[q],
                                                     density_0[q],
                                                     density_1[q]);

                  specific_heat =
                    calculate_point_property(filtered_phase_values[q],
                                             specific_heat_0[q],
                                             specific_heat_1[q]);
                }

              heat_in_domain +=
                (density * specific_heat * local_temperature_values[q]) *
                phase_coefficient * fe_values_ht.JxW(q);

            } // end loop on quadrature points
        }
    } // end loop on cells

  heat_in_domain = Utilities::MPI::sum(heat_in_domain, mpi_communicator);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Thermal energy in " << domain_name << " : "
                  << heat_in_domain << std::endl;
    }

  // Filling table
  this->heat_flux_table.add_value("thermal_energy_" + domain_name,
                                  heat_in_domain);
}

template void
HeatTransfer<2>::postprocess_thermal_energy_in_fluid<
  TrilinosWrappers::MPI::Vector>(
  const bool                           gather_vof,
  const Parameters::FluidIndicator     monitored_fluid,
  const std::string                    domain_name,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
HeatTransfer<3>::postprocess_thermal_energy_in_fluid<
  TrilinosWrappers::MPI::Vector>(
  const bool                           gather_vof,
  const Parameters::FluidIndicator     monitored_fluid,
  const std::string                    domain_name,
  const TrilinosWrappers::MPI::Vector &current_solution_fd);

template void
HeatTransfer<2>::postprocess_thermal_energy_in_fluid<
  TrilinosWrappers::MPI::BlockVector>(
  const bool                                gather_vof,
  const Parameters::FluidIndicator          monitored_fluid,
  const std::string                         domain_name,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template void
HeatTransfer<3>::postprocess_thermal_energy_in_fluid<
  TrilinosWrappers::MPI::BlockVector>(
  const bool                                gather_vof,
  const Parameters::FluidIndicator          monitored_fluid,
  const std::string                         domain_name,
  const TrilinosWrappers::MPI::BlockVector &current_solution_fd);

template <int dim>
void
HeatTransfer<dim>::write_heat_flux(const std::string domain_name)
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.heat_flux_output_name + "_" +
        domain_name + ".dat";
      std::ofstream output(filename.c_str());

      this->heat_flux_table.write_text(output);
    }
}

template <int dim>
std::pair<double, bool>
HeatTransfer<dim>::set_phase_coefficient(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const double                     phase_value_q)
{
  double phase_coefficient(0.);
  bool   point_is_in_postprocessed_fluid(false);

  switch (monitored_fluid)
    {
        default: {
          // Parameters::FluidIndicator::both
          phase_coefficient               = 1.;
          point_is_in_postprocessed_fluid = true;
          break;
        }
        case Parameters::FluidIndicator::fluid0: {
          if (gather_vof)
            {
              phase_coefficient = 1. - phase_value_q;
              if (phase_value_q < 0.5)
                point_is_in_postprocessed_fluid = true;
            }
          else
            {
              // In the case of monophase simulations, "fluid0" is
              // equivalent to "both" (all domain) calculation
              phase_coefficient               = 1.;
              point_is_in_postprocessed_fluid = true;
            }
          break;
        }
        case Parameters::FluidIndicator::fluid1: {
          if (gather_vof)
            {
              phase_coefficient = phase_value_q;
              if (phase_value_q > 0.5)
                point_is_in_postprocessed_fluid = true;
            }
          else
            {
              // Defensive only here, this has been checked already
              // in simulation_parameters.h
              throw std::logic_error(
                "Inconsistency in .prm!\n when VOF = false"
                "\n use (default value): set postprocessed fluid = both"
                "\n or: set postprocessed fluid = fluid 0");
            }
          break;
        }
    } // end switch on monitored fluid

  return std::make_pair(phase_coefficient, point_is_in_postprocessed_fluid);
}

template class HeatTransfer<2>;
template class HeatTransfer<3>;
