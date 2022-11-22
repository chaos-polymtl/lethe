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
  auto density_model =
    this->simulation_parameters.physical_properties_manager.get_density();
  auto specific_heat_model =
    this->simulation_parameters.physical_properties_manager.get_specific_heat();
  auto conductivity_model =
    this->simulation_parameters.physical_properties_manager
      .get_thermal_conductivity();
  std::map<field, double> field_values;

  double rho_cp = density_model->value(field_values) *
                  specific_heat_model->value(field_values);

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

              double h_cell = 0;
              if (dim == 2)
                h_cell =
                  std::sqrt(4. * cell->measure() / M_PI) /
                  this->simulation_parameters.fem_parameters.velocity_order;
              else if (dim == 3)
                h_cell =
                  pow(6 * cell->measure() / M_PI, 1. / 3.) /
                  this->simulation_parameters.fem_parameters.velocity_order;
              const double penalty_parameter =
                1. / std::pow(h_cell * h_cell, double(dim) / double(dim));
              const auto &dh_cell =
                typename DoFHandler<dim>::cell_iterator(*cell,
                                                        &this->dof_handler);
              dh_cell->get_dof_indices(heat_dof_indices);

              const auto pic = solid_ph->particles_in_cell(cell);
              Assert(pic.begin() == particle, ExcInternalError());
              for (const auto &p : pic)
                {
                  double      temperature = 0;
                  const auto &ref_q       = p.get_reference_location();
                  const auto &real_q      = p.get_location();
                  const auto &JxW         = p.get_properties()[0];

                  for (unsigned int k = 0; k < dofs_per_cell; ++k)
                    {
                      // Get the temperature at non-quadrature point (particle
                      // in fluid)
                      auto &evaluation_point = this->evaluation_point;
                      temperature += evaluation_point[heat_dof_indices[k]] *
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
                                rho_cp * penalty_parameter * beta *
                                this->fe->shape_value(i, ref_q) *
                                this->fe->shape_value(j, ref_q) * JxW;
                            }
                        }
                      else
                        {
                          // Regular residual
                          local_rhs(i) +=
                            -rho_cp * penalty_parameter * beta * temperature *
                              this->fe->shape_value(i, ref_q) * JxW +
                            rho_cp * penalty_parameter * beta *
                              solid_temperature->value(real_q, 0) *
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
HeatTransfer<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Laser heat source
  if (this->simulation_parameters.laser_parameters->activate_laser)
    {
      this->assemblers.push_back(
        std::make_shared<HeatTransferAssemblerLaser<dim>>(
          this->simulation_control,
          this->simulation_parameters.laser_parameters));
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
      scratch_data.enable_vof(dof_handler_vof->get_fe(),
                              *this->cell_quadrature,
                              *this->temperature_mapping);
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
            .use_time_average_velocity_field)
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
            .use_time_average_velocity_field)
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
                              *this->multiphysics->get_solution(PhysicsID::VOF),
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
      scratch_data.enable_vof(dof_handler_vof->get_fe(),
                              *this->cell_quadrature,
                              *this->temperature_mapping);
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
            .use_time_average_velocity_field)
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
            .use_time_average_velocity_field)
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
                              *this->multiphysics->get_solution(PhysicsID::VOF),
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

  std::string domain_name("");
  bool        postprocess_all_domain(false);

  if (monitored_fluid == Parameters::FluidIndicator::both ||
      (not(this->simulation_parameters.multiphysics.VOF) &&
       monitored_fluid == Parameters::FluidIndicator::fluid0))
    {
      domain_name            = "all_domain";
      postprocess_all_domain = true;
    }
  else
    {
      switch (monitored_fluid)
        {
          case Parameters::FluidIndicator::fluid0:
            {
              domain_name = "fluid_0";
              break;
            }
          case Parameters::FluidIndicator::fluid1:
            {
              domain_name = "fluid_1";
              break;
            }
          default:
            throw std::runtime_error("Unsupported number of fluids (>2)");
        }
    }

  // Temperature statistics
  if (simulation_parameters.post_processing.calculate_temperature_statistics)
    {
      if (postprocess_all_domain)
        {
          calculate_temperature_statistics_on_all_domain();
        }
      else
        {
          calculate_temperature_statistics_on_one_fluid(monitored_fluid,
                                                        domain_name);
        }

      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_temperature_statistics(domain_name);
    }

  // Heat flux calculation
  if (simulation_parameters.post_processing.calculate_heat_flux)
    {
      calculate_heat_flux_on_bc();

      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_heat_flux();
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
HeatTransfer<dim>::calculate_temperature_statistics_on_all_domain()
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  std::vector<double> local_temperature_values(n_q_points);

  FEValues<dim> fe_values_ht(*this->temperature_mapping,
                             *this->fe,
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  double volume_integral      = 0;
  double temperature_integral = 0;
  double minimum_temperature  = DBL_MAX;
  double maximum_temperature  = -DBL_MAX;

  // Calculate min, max and average
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              maximum_temperature =
                std::max(local_temperature_values[q], maximum_temperature);
              minimum_temperature =
                std::min(local_temperature_values[q], minimum_temperature);

              volume_integral += fe_values_ht.JxW(q);
              temperature_integral +=
                local_temperature_values[q] * fe_values_ht.JxW(q);
            }
        }
    }

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

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              temperature_variance_integral +=
                (local_temperature_values[q] - temperature_average) *
                (local_temperature_values[q] - temperature_average) *
                fe_values_ht.JxW(q);
            }
        }
    }

  temperature_variance_integral =
    Utilities::MPI::sum(temperature_variance_integral, mpi_communicator);
  double temperature_variance = temperature_variance_integral / volume_integral;
  double temperature_std_deviation = std::sqrt(temperature_variance);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Temperature statistics on all domain : " << std::endl;
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
HeatTransfer<dim>::calculate_temperature_statistics_on_one_fluid(
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  std::vector<double> local_temperature_values(n_q_points);

  FEValues<dim> fe_values_ht(*this->temperature_mapping,
                             *this->fe,
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Gather VOF information
  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics->get_dof_handler(PhysicsID::VOF);

  FEValues<dim> fe_values_vof(*this->temperature_mapping,
                              dof_handler_vof->get_fe(),
                              *this->cell_quadrature,
                              update_values | update_JxW_values);

  std::vector<double> phase_values(n_q_points);

  double phase_coefficient(0.);
  double volume_integral      = 0;
  double temperature_integral = 0;
  double minimum_temperature  = DBL_MAX;
  double maximum_temperature  = -DBL_MAX;

  // Calculate min, max and average
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          // Get VOF active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_vof(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_vof);

          // Gather VOF information
          fe_values_vof.reinit(cell_vof);
          fe_values_vof.get_function_values(
            *this->multiphysics->get_solution(PhysicsID::VOF), phase_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (monitored_fluid)
                {
                  case Parameters::FluidIndicator::fluid0:
                    {
                      phase_coefficient = 1. - phase_values[q];
                      if (phase_values[q] < 0.5)
                        {
                          maximum_temperature =
                            std::max(local_temperature_values[q],
                                     maximum_temperature);
                          minimum_temperature =
                            std::min(local_temperature_values[q],
                                     minimum_temperature);
                        }
                      break;
                    }
                  case Parameters::FluidIndicator::fluid1:
                    {
                      phase_coefficient = phase_values[q];
                      if (phase_values[q] > 0.5)
                        {
                          maximum_temperature =
                            std::max(local_temperature_values[q],
                                     maximum_temperature);
                          minimum_temperature =
                            std::min(local_temperature_values[q],
                                     minimum_temperature);
                        }
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
                } // end switch on monitored fluid

              volume_integral += fe_values_vof.JxW(q) * phase_coefficient;
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

          // Get VOF active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_vof(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_vof);

          // Gather VOF information
          fe_values_vof.reinit(cell_vof);
          fe_values_vof.get_function_values(
            *this->multiphysics->get_solution(PhysicsID::VOF), phase_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              switch (monitored_fluid)
                {
                  case Parameters::FluidIndicator::fluid0:
                    {
                      phase_coefficient = 1 - phase_values[q];
                      break;
                    }
                  case Parameters::FluidIndicator::fluid1:
                    {
                      phase_coefficient = phase_values[q];
                      break;
                    }
                  default:
                    throw std::runtime_error(
                      "Unsupported number of fluids (>2)");
                } // end switch on monitored fluid
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
        simulation_parameters.post_processing.temperature_output_name + "_" +
        domain_name + ".dat";
      std::ofstream output(filename.c_str());

      this->statistics_table.write_text(output);
    }
}

template <int dim>
void
HeatTransfer<dim>::calculate_heat_flux_on_bc()
{
  const FESystem<dim, dim> fe_ht = this->dof_handler.get_fe();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FESystem<dim, dim> fe_fd = dof_handler_fluid->get_fe();

  // Evaluate fluid properties
  auto density_model =
    this->simulation_parameters.physical_properties_manager.get_density();
  auto specific_heat_model =
    this->simulation_parameters.physical_properties_manager.get_specific_heat();
  auto conductivity_model =
    this->simulation_parameters.physical_properties_manager
      .get_thermal_conductivity();
  std::map<field, double> dummy_field_map;

  double rho_cp = density_model->value(dummy_field_map) *
                  specific_heat_model->value(dummy_field_map);

  double conductivity = conductivity_model->value(dummy_field_map);

  const unsigned int               n_q_points = this->face_quadrature->size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> velocity_values =
    std::vector<Tensor<1, dim>>(n_q_points);
  std::vector<Tensor<1, dim>> temperature_gradient =
    std::vector<Tensor<1, dim>>(n_q_points);
  std::vector<double> temperatures = std::vector<double>(n_q_points);

  Tensor<1, dim> normal_vector;
  double         heat_flux_bc;

  std::vector<double> heat_flux_vector(
    this->simulation_parameters.boundary_conditions.size);

  FEFaceValues<dim> fe_face_values_ht(*this->temperature_mapping,
                                      fe_ht,
                                      *this->face_quadrature,
                                      update_values | update_quadrature_points |
                                        update_gradients | update_JxW_values |
                                        update_normal_vectors);

  FEFaceValues<dim> fe_face_values_fd(*this->temperature_mapping,
                                      fe_fd,
                                      *this->face_quadrature,
                                      update_values | update_quadrature_points);

  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  TrilinosWrappers::MPI::Vector fluid_solution =
    *multiphysics->get_solution(PhysicsID::fluid_dynamics);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      unsigned int boundary_id =
        this->simulation_parameters.boundary_conditions.id[i_bc];
      heat_flux_bc = 0;

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              if (cell->at_boundary())
                {
                  for (const auto face : cell->face_indices())
                    {
                      if (cell->face(face)->at_boundary())
                        {
                          typename DoFHandler<dim>::active_cell_iterator
                            velocity_cell(&(*triangulation),
                                          cell->level(),
                                          cell->index(),
                                          dof_handler_fluid);
                          fe_face_values_ht.reinit(cell, face);
                          fe_face_values_fd.reinit(velocity_cell, face);

                          if (cell->face(face)->boundary_id() == boundary_id)
                            {
                              std::vector<Point<dim>> q_points =
                                fe_face_values_ht.get_quadrature_points();

                              fe_face_values_ht.get_function_gradients(
                                evaluation_point, temperature_gradient);

                              fe_face_values_ht.get_function_values(
                                evaluation_point, temperatures);

                              fe_face_values_fd[velocities].get_function_values(
                                fluid_solution, velocity_values);

                              for (unsigned int q = 0; q < n_q_points; q++)
                                {
                                  normal_vector =
                                    -fe_face_values_ht.normal_vector(q);

                                  heat_flux_bc +=
                                    (-conductivity * temperature_gradient[q] *
                                       normal_vector +
                                     temperatures[q] * rho_cp *
                                       velocity_values[q] * normal_vector) *
                                    fe_face_values_ht.JxW(q);
                                }
                            }
                        }
                    }
                }
            }
        }
      heat_flux_vector[i_bc] =
        Utilities::MPI::sum(heat_flux_bc, mpi_communicator);
    }

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Heat flux on heat transfer boundary conditions : "
                  << std::endl;
      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions.size;
           ++i_bc)
        this->pcout << "\t boundary " << i_bc << " : " << heat_flux_vector[i_bc]
                    << std::endl;
    }

  // Filling table
  this->heat_flux_table.add_value("time",
                                  this->simulation_control->get_current_time());
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    this->statistics_table.add_value("bc " + i_bc, heat_flux_vector[i_bc]);
}

template <int dim>
void
HeatTransfer<dim>::write_heat_flux()
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.post_processing.heat_flux_output_name + ".dat";
      std::ofstream output(filename.c_str());

      this->heat_flux_table.write_text(output);
    }
}

template class HeatTransfer<2>;
template class HeatTransfer<3>;
