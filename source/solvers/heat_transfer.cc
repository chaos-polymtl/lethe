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
#include <solvers/heat_transfer.h>


template <int dim>
void
HeatTransfer<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<true>(time_stepping_method);
}


template <int dim>
void
HeatTransfer<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<false>(time_stepping_method);
}


template <int dim>
template <bool assemble_matrix>
void
HeatTransfer<dim>::assemble_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  const double density = simulation_parameters.physical_properties.density;
  const double specific_heat =
    simulation_parameters.physical_properties.specific_heat;
  const double thermal_conductivity =
    simulation_parameters.physical_properties.thermal_conductivity;

  const double viscosity = simulation_parameters.physical_properties.viscosity;

  const double dynamic_viscosity = viscosity * density;

  const double rho_cp = density * specific_heat;

  auto &solution = present_solution;
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
        "SDIRK schemes are not supported by heat transfer physics");
    }


  auto &source_term = simulation_parameters.sourceTerm->heat_transfer_source;
  source_term.set_time(simulation_control->get_current_time());

  const QGauss<dim> quadrature_formula(fe.degree + 1);
  FEValues<dim>     fe_values_ht(fe,
                             quadrature_formula,
                             update_values | update_gradients |
                               update_quadrature_points | update_JxW_values);

  auto &evaluation_point = this->get_evaluation_point();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int                   n_q_points = quadrature_formula.size();
  std::vector<Tensor<1, dim>>          temperature_gradients(n_q_points);
  std::vector<double>                  source_term_values(n_q_points);


  const MappingQ<dim> mapping(
    fe.degree, simulation_parameters.fem_parameters.qmapping_all);

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  FEValues<dim> fe_values_flow(dof_handler_fluid->get_fe(),
                               quadrature_formula,
                               update_values | update_quadrature_points |
                                 update_gradients);

  // FaceValues for Robin boundary condition
  QGauss<dim - 1>   face_quadrature_formula(fe.degree + 1);
  FEFaceValues<dim> fe_face_values_ht(fe,
                                      face_quadrature_formula,
                                      update_values | update_quadrature_points |
                                        update_JxW_values);

  // Velocity values
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> velocity_gradient_values(n_q_points);

  std::vector<double> present_temperature_values(n_q_points);
  std::vector<double> present_face_temperature_values(
    face_quadrature_formula.size());

  // Values for backward Euler scheme
  std::vector<double> p1_temperature_values(n_q_points);
  std::vector<double> p2_temperature_values(n_q_points);
  std::vector<double> p3_temperature_values(n_q_points);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;
          fe_values_ht.reinit(cell);

          fe_values_ht.get_function_gradients(evaluation_point,
                                              temperature_gradients);


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
          fe_values_ht.get_function_values(solution,
                                           present_temperature_values);

          // Gather the previous time steps for heat transfer depending on
          // the number of stages of the time integration method
          if (time_stepping_method !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            fe_values_ht.get_function_values(this->solution_m1,
                                             p1_temperature_values);

          if (time_stepping_method_has_two_stages(time_stepping_method))
            fe_values_ht.get_function_values(this->solution_m2,
                                             p2_temperature_values);

          if (time_stepping_method_has_three_stages(time_stepping_method))
            fe_values_ht.get_function_values(this->solution_m3,
                                             p3_temperature_values);

          source_term.value_list(fe_values_ht.get_quadrature_points(),
                                 source_term_values);


          // assembling local matrix and right hand side
          for (const unsigned int q : fe_values_ht.quadrature_point_indices())
            {
              for (const unsigned int i : fe_values_ht.dof_indices())
                {
                  if (assemble_matrix)
                    {
                      for (const unsigned int j : fe_values_ht.dof_indices())
                        {
                          // Weak form for : - k * laplacian T + rho * cp * u *
                          //                      gradT - f - grad(u)*grad(u) =0
                          cell_matrix(i, j) +=
                            (thermal_conductivity *
                               fe_values_ht.shape_grad(i, q) *
                               fe_values_ht.shape_grad(j, q) +
                             rho_cp * fe_values_ht.shape_value(i, q) *
                               velocity_values[q] *
                               fe_values_ht.shape_grad(j,
                                                       q)) *
                            fe_values_ht.JxW(q); // JxW

                          // Mass matrix for transient simulation
                          if (is_bdf(time_stepping_method))
                            cell_matrix(i, j) +=
                              rho_cp * fe_values_ht.shape_value(j, q) *
                              fe_values_ht.shape_value(i, q) * bdf_coefs[0] *
                              fe_values_ht.JxW(q);
                        }
                    }

                  // rhs for : - k * laplacian T + rho * cp * u * grad T - f
                  // -grad(u)*grad(u) = 0
                  cell_rhs(i) -=
                    (thermal_conductivity * fe_values_ht.shape_grad(i, q) *
                       temperature_gradients[q] +
                     +density * specific_heat * fe_values_ht.shape_value(i, q) *
                       velocity_values[q] * temperature_gradients[q] -
                     source_term_values[q] * fe_values_ht.shape_value(i, q) -
                     dynamic_viscosity * fe_values_ht.shape_value(i, q) *
                       scalar_product(velocity_gradient_values[q] +
                                        transpose(velocity_gradient_values[q]),
                                      transpose(velocity_gradient_values[q]))) *
                    fe_values_ht.JxW(q); // JxW

                  // Residual associated with BDF schemes
                  if (time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::bdf1 ||
                      time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::steady_bdf)
                    cell_rhs(i) -=
                      rho_cp *
                      (bdf_coefs[0] * present_temperature_values[q] +
                       bdf_coefs[1] * p1_temperature_values[q]) *
                      fe_values_ht.shape_value(i, q) *
                      fe_values_ht.JxW(q); // *phi_u[i]*JxW

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    cell_rhs(i) -=
                      rho_cp *
                      (bdf_coefs[0] * present_temperature_values[q] +
                       bdf_coefs[1] * p1_temperature_values[q] +
                       bdf_coefs[2] * p2_temperature_values[q]) *
                      fe_values_ht.shape_value(i, q) *
                      fe_values_ht.JxW(q); // *phi_u[i]*JxW

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    cell_rhs(i) -=
                      rho_cp *
                      (bdf_coefs[0] * present_temperature_values[q] +
                       bdf_coefs[1] * p1_temperature_values[q] +
                       bdf_coefs[2] * p2_temperature_values[q] +
                       bdf_coefs[3] * p3_temperature_values[q]) *
                      fe_values_ht.shape_value(i, q) *
                      fe_values_ht.JxW(q); // *phi_u[i]*JxW
                }

            } // end loop on quadrature points

          // Robin boundary condition, loop on faces (Newton's cooling law)
          // implementation similar to deal.ii step-7
          for (unsigned int i_bc = 0;
               i_bc < simulation_parameters.boundary_conditions_ht.size;
               ++i_bc)
            {
              if (this->simulation_parameters.boundary_conditions_ht
                    .type[i_bc] == BoundaryConditions::BoundaryType::convection)
                {
                  const double h =
                    simulation_parameters.boundary_conditions_ht.h[i_bc];
                  const double T_inf =
                    simulation_parameters.boundary_conditions_ht.Tinf[i_bc];
                  if (cell->is_locally_owned())
                    {
                      for (unsigned int face = 0;
                           face < GeometryInfo<dim>::faces_per_cell;
                           face++)
                        {
                          if (cell->face(face)->at_boundary() &&
                              (cell->face(face)->boundary_id() ==
                               simulation_parameters.boundary_conditions_ht
                                 .id[i_bc]))
                            {
                              fe_face_values_ht.reinit(cell, face);
                              fe_face_values_ht.get_function_values(
                                solution, present_face_temperature_values);
                              {
                                for (const unsigned int q :
                                     fe_face_values_ht
                                       .quadrature_point_indices())
                                  {
                                    const double JxW = fe_face_values_ht.JxW(q);
                                    for (const unsigned int i :
                                         fe_values_ht.dof_indices())
                                      {
                                        if (assemble_matrix)
                                          {
                                            for (const unsigned int j :
                                                 fe_values_ht.dof_indices())
                                              {
                                                // Weak form modification
                                                cell_matrix(i, j) +=
                                                  fe_face_values_ht.shape_value(
                                                    i, q) *
                                                  fe_face_values_ht.shape_value(
                                                    j, q) *
                                                  h * JxW;
                                              }
                                          }
                                        // Residual
                                        cell_rhs(i) -=
                                          fe_face_values_ht.shape_value(i, q) *
                                          h *
                                          (present_face_temperature_values[q] -
                                           T_inf) *
                                          JxW;
                                      }
                                  }
                              }
                            }
                        }
                    }
                }
            } // end loop for Robin condition

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
HeatTransfer<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  data_out.add_data_vector(dof_handler, present_solution, "temperature");
}

template <int dim>
double
HeatTransfer<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();


  QGauss<dim>         quadrature_formula(fe.degree + 2);
  const MappingQ<dim> mapping(
    fe.degree, simulation_parameters.fem_parameters.qmapping_all);
  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);



  const unsigned int dofs_per_cell =
    fe.dofs_per_cell; // This gives you dofs per cell

  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

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
HeatTransfer<dim>::finish_simulation()
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
      error_table.set_scientific("error_temperature", true);
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
HeatTransfer<dim>::percolate_time_vectors()
{
  solution_m3 = solution_m2;
  solution_m2 = solution_m1;
  solution_m1 = present_solution;
}

template <int dim>
void
HeatTransfer<dim>::finish_time_step()
{
  percolate_time_vectors();
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
}

template <int dim>
void
HeatTransfer<dim>::pre_mesh_adaptation()
{
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(solution_m3);
}

template <int dim>
void
HeatTransfer<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_communicator();


  // Set up the vectors for the transfer
  TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);
  TrilinosWrappers::MPI::Vector tmp_m1(locally_owned_dofs, mpi_communicator);
  TrilinosWrappers::MPI::Vector tmp_m2(locally_owned_dofs, mpi_communicator);
  TrilinosWrappers::MPI::Vector tmp_m3(locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);
  solution_transfer_m1.interpolate(tmp_m1);
  solution_transfer_m2.interpolate(tmp_m2);
  solution_transfer_m3.interpolate(tmp_m3);

  // Distribute constraints
  nonzero_constraints.distribute(tmp);
  nonzero_constraints.distribute(tmp_m1);
  nonzero_constraints.distribute(tmp_m2);
  nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  present_solution = tmp;
  solution_m1      = tmp_m1;
  solution_m2      = tmp_m2;
  solution_m3      = tmp_m3;
}



template <int dim>
void
HeatTransfer<dim>::setup_dofs()
{
  dof_handler.distribute_dofs(fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  auto mpi_communicator = triangulation->get_communicator();


  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  present_solution.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);

  // Previous solutions for transient schemes
  solution_m1.reinit(locally_owned_dofs,
                     locally_relevant_dofs,
                     mpi_communicator);
  solution_m2.reinit(locally_owned_dofs,
                     locally_relevant_dofs,
                     mpi_communicator);
  solution_m3.reinit(locally_owned_dofs,
                     locally_relevant_dofs,
                     mpi_communicator);

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
}

template <int dim>
void
HeatTransfer<dim>::set_initial_conditions()
{
  MappingQ<dim> mapping(fe.degree);
  VectorTools::interpolate(mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->temperature,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  finish_time_step();
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



template class HeatTransfer<2>;
template class HeatTransfer<3>;
