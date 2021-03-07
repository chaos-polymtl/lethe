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
#include <solvers/free_surface.h>


template <int dim>
void
FreeSurface<dim>::assemble_matrix_and_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<true>(time_stepping_method);
}


template <int dim>
void
FreeSurface<dim>::assemble_rhs(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
  assemble_system<false>(time_stepping_method);
}


template <int dim>
template <bool assemble_matrix>
void
FreeSurface<dim>::assemble_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method)
{
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
  const double dt   = time_steps_vector[0];
  const double sdt2 = std::pow(1. / dt, 2);

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


  // Free surface FEValues initialization
  // TODO clean up updated values after full implementation
  FEValues<dim> fe_values_fs(*fe,
                             *this->cell_quadrature,
                             update_values | update_gradients |
                               update_quadrature_points | update_JxW_values |
                               update_hessians);

  auto &evaluation_point = this->get_evaluation_point();

  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const unsigned int  n_q_points = this->cell_quadrature->size();
  std::vector<double> source_term_values(n_q_points);

  // Fluid FEValues information gathering
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  FEValues<dim> fe_values_flow(dof_handler_fluid->get_fe(),
                               *this->cell_quadrature,
                               update_values | update_quadrature_points |
                                 update_gradients);

  // Shape functions and gradients
  std::vector<double>         phi_phase(dofs_per_cell);
  std::vector<Tensor<1, dim>> grad_phi_phase(dofs_per_cell);
  std::vector<Tensor<2, dim>> hess_phi_phase(dofs_per_cell);
  std::vector<double>         laplacian_phi_phase(dofs_per_cell);


  // Velocity values
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> velocity_gradient_values(n_q_points);

  std::vector<double>         present_phase_values(n_q_points);
  std::vector<Tensor<1, dim>> phase_gradients(n_q_points);
  std::vector<double>         present_phase_laplacians(n_q_points);

  // Values for backward Euler scheme
  std::vector<double> p1_phase_values(n_q_points);
  std::vector<double> p2_phase_values(n_q_points);
  std::vector<double> p3_phase_values(n_q_points);

  // Element size
  double h;

  double velocity_fem_degree =
    simulation_parameters.fem_parameters.velocity_order;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_matrix = 0;
          cell_rhs    = 0;

          fe_values_fs.reinit(cell);

          fe_values_fs.get_function_gradients(evaluation_point,
                                              phase_gradients);

          // Element size (NB : dim is implicitly converted to double)
          h = pow(2. * dim * cell->measure() / M_PI, 1. / dim) /
              velocity_fem_degree;

          // Gather flow values
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

          // Gather free surface present value
          fe_values_fs.get_function_values(evaluation_point,
                                           present_phase_values);

          // Gather present laplacian
          fe_values_fs.get_function_laplacians(evaluation_point,
                                               present_phase_laplacians);

          // Gather the previous time steps for heat transfer depending on
          // the number of stages of the time integration method
          if (time_stepping_method !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_fs.get_function_values(this->solution_m1,
                                               p1_phase_values);
              // TODO see if needed
              //              fe_values_fs.get_function_gradients(this->solution_m1,
              //                                                  p1_phase_gradients);
            }

          if (time_stepping_method_has_two_stages(time_stepping_method))
            {
              fe_values_fs.get_function_values(this->solution_m2,
                                               p2_phase_values);

              //              fe_values_fs.get_function_gradients(this->solution_m2,
              //                                                  p2_phase_gradients);
            }

          if (time_stepping_method_has_three_stages(time_stepping_method))
            {
              fe_values_fs.get_function_values(this->solution_m3,
                                               p3_phase_values);

              //              fe_values_fs.get_function_gradients(this->solution_m3,
              //                                                  p3_phase_gradients);
            }


          // assembling local matrix and right hand side
          for (const unsigned int q : fe_values_fs.quadrature_point_indices())
            {
              // Store JxW in local variable for faster access
              const double JxW = fe_values_fs.JxW(q);

              const auto velocity      = velocity_values[q];
              const auto present_phase = present_phase_values[q];

              // Calculation of the magnitude of the velocity for the
              // stabilization parameter and the compression term for the phase
              // indicator
              const double u_mag = std::max(velocity.norm(), 1e-12);

              // Compression ratio //TODO use for compression term
              //              const double cr =
              //                std::max(velocity.norm(), 1e-12) /
              //                present_phase_values[q];

              // Calculation of the GLS stabilization parameter. The
              // stabilization parameter used is different if the simulation is
              // steady or unsteady. In the unsteady case it includes the value
              // of the time-step. Hypothesis : advection dominated problem
              // (Pe>3) [Bochev et al., Stability of the SUPG finite element
              // method for transient advection-diffusion problems, CMAME 2004]
              const double tau =
                is_steady(time_stepping_method) ?
                  h / (2. * u_mag) :
                  1. / std::sqrt(std::pow(2. * u_mag / h, 2) + sdt2);

              // Gather the shape functions and their gradient
              for (unsigned int k : fe_values_fs.dof_indices())
                {
                  phi_phase[k]      = fe_values_fs.shape_value(k, q);
                  grad_phi_phase[k] = fe_values_fs.shape_grad(k, q);
                  hess_phi_phase[k] = fe_values_fs.shape_hessian(k, q);

                  laplacian_phi_phase[k] = trace(hess_phi_phase[k]);
                }

              for (const unsigned int i : fe_values_fs.dof_indices())
                {
                  const auto phi_phase_i      = phi_phase[i];
                  const auto grad_phi_phase_i = grad_phi_phase[i];


                  if (assemble_matrix)
                    {
                      for (const unsigned int j : fe_values_fs.dof_indices())
                        {
                          const auto phi_phase_j      = phi_phase[j];
                          const auto grad_phi_phase_j = grad_phi_phase[j];
                          //                          const auto
                          //                          laplacian_phi_phase_j =
                          //                            laplacian_phi_phase[j];



                          // Weak form for advection: u * nabla(phase) = 0
                          // TODO add compression term
                          cell_matrix(i, j) +=
                            (phi_phase_i * velocity * grad_phi_phase_j) * JxW;

                          // Weak form for compression term:
                          // div(phase(1-phase)*nabla(phase)*cr)

                          // Strong Jacobian associated with the GLS
                          // stabilization
                          auto strong_jacobian = velocity * grad_phi_phase_j;

                          // Mass matrix for transient simulation
                          if (is_bdf(time_stepping_method))
                            {
                              cell_matrix(i, j) +=
                                phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW;

                              strong_jacobian += phi_phase_j * bdf_coefs[0];
                            }

                          // Addition to the cell matrix for GLS stabilization
                          cell_matrix(i, j) +=
                            tau * strong_jacobian *
                            (grad_phi_phase_i * velocity_values[q]) * JxW;
                        }
                    }

                  // rhs for : u * nabla(phase) = 0
                  // TODO add compression term
                  cell_rhs(i) -=
                    (phi_phase_i * velocity_values[q] * phase_gradients[q]) *
                    JxW;

                  // Calculate the strong residual for GLS stabilization
                  auto strong_residual =
                    velocity_values[q] * phase_gradients[q];

                  // Residual associated with BDF schemes
                  if (time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::bdf1 ||
                      time_stepping_method == Parameters::SimulationControl::
                                                TimeSteppingMethod::steady_bdf)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_phase +
                                      bdf_coefs[1] * p1_phase_values[q]) *
                                     phi_phase_i * JxW;

                      strong_residual += bdf_coefs[0] * present_phase +
                                         bdf_coefs[1] * p1_phase_values[q];
                    }

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_phase +
                                      bdf_coefs[1] * p1_phase_values[q] +
                                      bdf_coefs[2] * p2_phase_values[q]) *
                                     phi_phase_i * JxW;

                      strong_residual += bdf_coefs[0] * present_phase +
                                         bdf_coefs[1] * p1_phase_values[q] +
                                         bdf_coefs[2] * p2_phase_values[q];
                    }

                  if (time_stepping_method ==
                      Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                    {
                      cell_rhs(i) -= (bdf_coefs[0] * present_phase +
                                      bdf_coefs[1] * p1_phase_values[q] +
                                      bdf_coefs[2] * p2_phase_values[q] +
                                      bdf_coefs[3] * p3_phase_values[q]) *
                                     phi_phase_i * JxW;

                      strong_residual += (bdf_coefs[0] * present_phase +
                                          bdf_coefs[1] * p1_phase_values[q] +
                                          bdf_coefs[2] * p2_phase_values[q] +
                                          bdf_coefs[3] * p3_phase_values[q]);
                    }


                  cell_rhs(i) -= tau *
                                 (strong_residual *
                                  (grad_phi_phase_i * velocity_values[q])) *
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
FreeSurface<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  data_out.add_data_vector(dof_handler, present_solution, "phase");
}

template <int dim>
double
FreeSurface<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values_fs(*this->mapping,
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
void
FreeSurface<dim>::finish_simulation()
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
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
FreeSurface<dim>::percolate_time_vectors()
{
  solution_m3 = solution_m2;
  solution_m2 = solution_m1;
  solution_m1 = present_solution;
}

template <int dim>
void
FreeSurface<dim>::finish_time_step()
{
  percolate_time_vectors();
}

template <int dim>
void
FreeSurface<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
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
}

template <int dim>
void
FreeSurface<dim>::pre_mesh_adaptation()
{
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(solution_m3);
}

template <int dim>
void
FreeSurface<dim>::post_mesh_adaptation()
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
FreeSurface<dim>::write_checkpoint()
{
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;

  sol_set_transfer.push_back(&present_solution);
  sol_set_transfer.push_back(&solution_m1);
  sol_set_transfer.push_back(&solution_m2);
  sol_set_transfer.push_back(&solution_m3);
  solution_transfer.prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
FreeSurface<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading free surface checkpoint" << std::endl;

  std::vector<TrilinosWrappers::MPI::Vector *> input_vectors(4);
  TrilinosWrappers::MPI::Vector distributed_system(locally_owned_dofs,
                                                   mpi_communicator);
  TrilinosWrappers::MPI::Vector distributed_system_m1(locally_owned_dofs,
                                                      mpi_communicator);
  TrilinosWrappers::MPI::Vector distributed_system_m2(locally_owned_dofs,
                                                      mpi_communicator);
  TrilinosWrappers::MPI::Vector distributed_system_m3(locally_owned_dofs,
                                                      mpi_communicator);

  input_vectors[0] = &distributed_system;
  input_vectors[1] = &distributed_system_m1;
  input_vectors[2] = &distributed_system_m2;
  input_vectors[3] = &distributed_system_m3;

  solution_transfer.deserialize(input_vectors);

  present_solution = distributed_system;
  solution_m1      = distributed_system_m1;
  solution_m2      = distributed_system_m2;
  solution_m3      = distributed_system_m3;
}


template <int dim>
void
FreeSurface<dim>::setup_dofs()
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
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of free surface degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the free surface dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::free_surface, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::free_surface, &this->present_solution);
  // the fluid at present iteration is solved before the free surface, and after
  // percolate is called for the previous iteration, so at the time the getter
  // is called solution_m2 = solution_m1
  // TODO deactivated for now (inertia is considered with a constant density),
  // see if needed / to be debugged
  multiphysics->set_solution_m1(PhysicsID::free_surface, &this->solution_m1);
}

template <int dim>
void
FreeSurface<dim>::set_initial_conditions()
{
  VectorTools::interpolate(
    *this->mapping,
    dof_handler,
    simulation_parameters.initial_condition->free_surface,
    newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;

  finish_time_step();
}

template <int dim>
void
FreeSurface<dim>::solve_linear_system(const bool initial_step,
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
      this->pcout << "  Free Surface : " << std::endl
                  << "  -Tolerance of iterative solver is : "
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
      this->pcout << "  Free Surface : " << std::endl
                  << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  // Update constraints and newton vectors
  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template class FreeSurface<2>;
template class FreeSurface<3>;
