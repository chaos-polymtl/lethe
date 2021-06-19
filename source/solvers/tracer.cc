#include <core/bdf.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/tracer.h>

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

#include <deal.II/numerics/vector_tools.h>


template <int dim>
void
Tracer<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<true>(time_stepping_method);
}


template <int dim>
void
Tracer<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<false>(time_stepping_method);
}


template <int dim>
template <bool assemble_matrix>
void
Tracer<dim>::assemble_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  const double tracer_diffusivity =
    simulation_parameters.physical_properties.tracer_diffusivity;

  if (assemble_matrix)
    system_matrix = 0;
  system_rhs = 0;

  // Vector for the BDF coefficients
  // The coefficients are stored in the following fashion :
  // 0 - n+1
  // 1 - n
  // 2 - n-1
  // 3 - n-2
  std::vector<double> time_steps_vector =
    simulation_control->get_time_steps_vector();

  // Time steps and inverse time steps which is used for numerous calculations
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;


  Vector<double> bdf_coefs;

  if (time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
      time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    bdf_coefs = bdf_coefficients(1, time_steps_vector);

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    bdf_coefs = bdf_coefficients(2, time_steps_vector);

  if (time_stepping_method ==
      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    bdf_coefs = bdf_coefficients(3, time_steps_vector);

  if (time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1 ||
      time_stepping_method ==
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1)
    {
      throw std::runtime_error(
        "SDIRK schemes are not supported by tracer physics");
    }


  auto &source_term = simulation_parameters.sourceTerm->tracer_source;
  source_term.set_time(simulation_control->get_current_time());

  FEValues<dim> fe_values_tracer(*mapping,
                                 *fe,
                                 *cell_quadrature,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values | update_hessians);

  auto &evaluation_point = this->get_evaluation_point();

  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int                   n_q_points = cell_quadrature->size();
  std::vector<double>                  source_term_values(n_q_points);


  const MappingQ<dim> mapping(
    fe->degree, simulation_parameters.fem_parameters.qmapping_all);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  FEValues<dim> fe_values_flow(dof_handler_fluid->get_fe(),
                               *cell_quadrature,
                               update_values | update_quadrature_points |
                                 update_gradients);

  // Shape functions and gradients
  std::vector<double>         phi_T(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell);
  std::vector<Tensor<2, dim>> hess_phi_T(dofs_per_cell);
  std::vector<double>         laplacian_phi_T(dofs_per_cell);


  // Velocity values
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> velocity_gradient_values(n_q_points);

  std::vector<double>         present_tracer_values(n_q_points);
  std::vector<Tensor<1, dim>> tracer_gradients(n_q_points);
  std::vector<double>         present_tracer_laplacians(n_q_points);

  // Values for backward Euler scheme
  std::vector<double> p1_tracer_values(n_q_points);
  std::vector<double> p2_tracer_values(n_q_points);
  std::vector<double> p3_tracer_values(n_q_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;
          double h    = 0;

          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) / fe->degree;
          else if (dim == 3)
            h = pow(6 * cell->measure() / M_PI, 1. / 3.) / fe->degree;

          fe_values_tracer.reinit(cell);

          fe_values_tracer.get_function_gradients(evaluation_point,
                                                  tracer_gradients);


          typename DoFHandler<dim>::active_cell_iterator velocity_cell(
            &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

          fe_values_flow.reinit(velocity_cell);

          if (multiphysics->fluid_dynamics_is_block())
            {
              fe_values_flow[velocities].get_function_values(
                *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
                velocity_values);
              fe_values_flow[velocities].get_function_gradients(
                *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
                velocity_gradient_values);
            }
          else
            {
              fe_values_flow[velocities].get_function_values(
                *multiphysics->get_solution(PhysicsID::fluid_dynamics),
                velocity_values);
              fe_values_flow[velocities].get_function_gradients(
                *multiphysics->get_solution(PhysicsID::fluid_dynamics),
                velocity_gradient_values);
            }

          // Gather present value
          fe_values_tracer.get_function_values(evaluation_point,
                                               present_tracer_values);

          // Gather present laplacian
          fe_values_tracer.get_function_laplacians(evaluation_point,
                                                   present_tracer_laplacians);

          // Gather the previous time steps for tracer depending on
          // the number of stages of the time integration method
          if (time_stepping_method !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_tracer.get_function_values(previous_solutions[0],
                                                   p1_tracer_values);
            }

          if (time_stepping_method_uses_two_previous_solutions(
                time_stepping_method))
            {
              fe_values_tracer.get_function_values(previous_solutions[1],
                                                   p2_tracer_values);
            }

          if (time_stepping_method_uses_three_previous_solutions(
                time_stepping_method))
            {
              fe_values_tracer.get_function_values(previous_solutions[2],
                                                   p3_tracer_values);
            }

          source_term.value_list(fe_values_tracer.get_quadrature_points(),
                                 source_term_values);


          // assembling local matrix and right hand side
          for (const unsigned int q :
               fe_values_tracer.quadrature_point_indices())
            {
              // Store JxW in local variable for faster access
              const double JxW = fe_values_tracer.JxW(q);

              const auto velocity = velocity_values[q];


              // Calculation of the magnitude of the velocity for the
              // stabilization parameter
              const double u_mag = std::max(velocity.norm(), 1e-12);

              // Calculation of the GLS stabilization parameter. The
              // stabilization parameter used is different if the simulation is
              // steady or unsteady. In the unsteady case it includes the value
              // of the time-step
              const double tau =
                is_steady(time_stepping_method) ?
                  1. / std::sqrt(
                         std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * tracer_diffusivity / (h * h), 2)) :
                  1. / std::sqrt(
                         std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * tracer_diffusivity / (h * h), 2));

              // Gather the shape functions and their gradient
              for (unsigned int k : fe_values_tracer.dof_indices())
                {
                  phi_T[k]      = fe_values_tracer.shape_value(k, q);
                  grad_phi_T[k] = fe_values_tracer.shape_grad(k, q);
                  hess_phi_T[k] = fe_values_tracer.shape_hessian(k, q);

                  laplacian_phi_T[k] = trace(hess_phi_T[k]);
                }

              for (const unsigned int i : fe_values_tracer.dof_indices())
                {
                  const auto phi_T_i      = phi_T[i];
                  const auto grad_phi_T_i = grad_phi_T[i];


                  if (assemble_matrix)
                    {
                      for (const unsigned int j :
                           fe_values_tracer.dof_indices())
                        {
                          const auto phi_T_j           = phi_T[j];
                          const auto grad_phi_T_j      = grad_phi_T[j];
                          const auto laplacian_phi_T_j = laplacian_phi_T[j];


                          // Weak form : - D * laplacian T +  u * gradT - f=0
                          cell_matrix(i, j) +=
                            (tracer_diffusivity * grad_phi_T_i * grad_phi_T_j +
                             phi_T_i * velocity * grad_phi_T_j) *
                            JxW;

                          auto strong_jacobian =
                            velocity * grad_phi_T_j -
                            tracer_diffusivity * laplacian_phi_T_j;

                          // Mass matrix for transient simulation
                          if (is_bdf(time_stepping_method))
                            {
                              cell_matrix(i, j) +=
                                phi_T_j * phi_T_i * bdf_coefs[0] * JxW;

                              strong_jacobian += phi_T_j * bdf_coefs[0];
                            }

                          cell_matrix(i, j) +=
                            tau * strong_jacobian *
                            (grad_phi_T_i * velocity_values[q]) * JxW;
                        }
                    }

                  // rhs for : - D * laplacian T +  u * grad T - f=0
                  cell_rhs(i) -=
                    (tracer_diffusivity * grad_phi_T_i * tracer_gradients[q] +
                     phi_T_i * velocity_values[q] * tracer_gradients[q] -
                     source_term_values[q] * phi_T_i) *
                    JxW;

                  // Calculate the strong residual for GLS stabilization
                  auto strong_residual =
                    velocity_values[q] * tracer_gradients[q] -
                    tracer_diffusivity * present_tracer_laplacians[q];



                  // Residual associated with BDF schemes
                  if (time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::bdf1 ||
                      time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::steady_bdf)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_tracer_values[q] +
                                      bdf_coefs[1] * p1_tracer_values[q]) *
                                     phi_T_i * JxW;

                      strong_residual +=
                        (bdf_coefs[0] * present_tracer_values[q] +
                         bdf_coefs[1] * p1_tracer_values[q]);
                    }

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_tracer_values[q] +
                                      bdf_coefs[1] * p1_tracer_values[q] +
                                      bdf_coefs[2] * p2_tracer_values[q]) *
                                     phi_T_i * JxW;

                      strong_residual +=
                        (bdf_coefs[0] * present_tracer_values[q] +
                         bdf_coefs[1] * p1_tracer_values[q] +
                         bdf_coefs[2] * p2_tracer_values[q]);
                    }

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_tracer_values[q] +
                                      bdf_coefs[1] * p1_tracer_values[q] +
                                      bdf_coefs[2] * p2_tracer_values[q] +
                                      bdf_coefs[3] * p3_tracer_values[q]) *
                                     phi_T_i * JxW;

                      strong_residual +=
                        (bdf_coefs[0] * present_tracer_values[q] +
                         bdf_coefs[1] * p1_tracer_values[q] +
                         bdf_coefs[2] * p2_tracer_values[q] +
                         bdf_coefs[3] * p3_tracer_values[q]);
                    }


                  cell_rhs(i) -=
                    tau *
                    (strong_residual * (grad_phi_T_i * velocity_values[q])) *
                    JxW;
                }

            } // end loop on quadrature points

          // transfer cell contribution into global objects
          cell->get_dof_indices(local_dof_indices);
          zero_constraints.distribute_local_to_global(cell_matrix,
                                                      cell_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
        } // end loop active cell
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}


template <int dim>
void
Tracer<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
#if (DEAL_II_VERSION_MAJOR < 10)
  data_out.add_data_vector(dof_handler, present_solution, "tracer");
#else
  local_evaluation_point = present_solution;
  data_out.add_data_vector(dof_handler, local_evaluation_point, "tracer");
#endif
}

template <int dim>
double
Tracer<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell =
    fe->dofs_per_cell; // This gives you dofs per cell

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = cell_quadrature->size();

  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->tracer;
  exact_solution.set_time(simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_scalar_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

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
Tracer<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_communicator();
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
      error_table.set_scientific("error_tracer", true);
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
Tracer<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
Tracer<dim>::finish_time_step()
{
  percolate_time_vectors();
}

template <int dim>
void
Tracer<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      double tracer_error = calculate_L2_error();

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());
      error_table.add_value("error_tracer", tracer_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error tracer : " << tracer_error << std::endl;
        }
    }

  if (simulation_parameters.post_processing.calculate_tracer_statistics)
    {
      calculate_tracer_statistics();
      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_tracer_statistics();
    }
}

template <int dim>
void
Tracer<dim>::calculate_tracer_statistics()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell =
    fe->dofs_per_cell; // This gives you dofs per cell

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int  n_q_points = cell_quadrature->size();
  std::vector<double> q_tracer_values(n_q_points);

  double volume_integral  = 0;
  double max_tracer_value = DBL_MIN;
  double min_tracer_value = DBL_MAX;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_tracer_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              volume_integral += q_tracer_values[q] * fe_values.JxW(q);
              max_tracer_value = std::max(q_tracer_values[q], max_tracer_value);
              min_tracer_value = std::min(q_tracer_values[q], min_tracer_value);
            }
        }
    }
  volume_integral      = Utilities::MPI::sum(volume_integral, mpi_communicator);
  double global_volume = GridTools::volume(*triangulation, *mapping);
  double tracer_average = volume_integral / global_volume;

  double variance_integral = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_tracer_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              variance_integral += (q_tracer_values[q] - tracer_average) *
                                   (q_tracer_values[q] - tracer_average) *
                                   fe_values.JxW(q);
            }
        }
    }

  variance_integral = Utilities::MPI::sum(variance_integral, mpi_communicator);
  double tracer_variance      = variance_integral / global_volume;
  double tracer_std_deviation = std::sqrt(tracer_variance);

  this->pcout << "Tracer statistics : " << std::endl;
  this->pcout << "\t     Min : " << min_tracer_value << std::endl;
  this->pcout << "\t     Max : " << max_tracer_value << std::endl;
  this->pcout << "\t Average : " << tracer_average << std::endl;
  this->pcout << "\t Std-Dev : " << tracer_std_deviation << std::endl;

  statistics_table.add_value("time", simulation_control->get_current_time());
  statistics_table.add_value("min", min_tracer_value);
  statistics_table.add_value("max", max_tracer_value);
  statistics_table.add_value("average", tracer_average);
  statistics_table.add_value("std-dev", tracer_std_deviation);
}

template <int dim>
void
Tracer<dim>::write_tracer_statistics()
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.post_processing.tracer_output_name + ".dat";
      std::ofstream output(filename.c_str());

      statistics_table.write_text(output);
    }
}

template <int dim>
void
Tracer<dim>::pre_mesh_adaptation()
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
Tracer<dim>::post_mesh_adaptation()
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
Tracer<dim>::write_checkpoint()
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
Tracer<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading tracer checkpoint" << std::endl;

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
Tracer<dim>::setup_dofs()
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
         i_bc < this->simulation_parameters.boundary_conditions_tracer.size;
         ++i_bc)
      {
        // Dirichlet condition : imposed temperature at i_bc
        if (this->simulation_parameters.boundary_conditions_tracer.type[i_bc] ==
            BoundaryConditions::BoundaryType::tracer_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_tracer.id[i_bc],
              *this->simulation_parameters.boundary_conditions_tracer
                 .tracer[i_bc],
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
         i_bc < this->simulation_parameters.boundary_conditions_tracer.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions_tracer.type[i_bc] ==
            BoundaryConditions::BoundaryType::tracer_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_tracer.id[i_bc],
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

  this->pcout << "   Number of tracer degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the tracer dof_handler and present solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::tracer, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::tracer, &this->present_solution);
}

template <int dim>
void
Tracer<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->tracer,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  finish_time_step();
}

template <int dim>
void
Tracer<dim>::solve_linear_system(const bool initial_step,
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



template class Tracer<2>;
template class Tracer<3>;
