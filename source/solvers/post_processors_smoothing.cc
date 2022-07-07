#include <solvers/post_processors_smoothing.h>

template <int dim, typename VectorType>
PostProcessorSmoothing<dim, VectorType>::PostProcessorSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points,
  const MPI_Comm &          mpi_communicator)
  : fe_q(1)
  , dof_handler(*triangulation)
  , simulation_parameters(simulation_parameters)
  , number_quadrature_points(number_quadrature_points)
  , mpi_communicator(mpi_communicator)
{
  system_matrix = std::make_shared<TrilinosWrappers::SparseMatrix>();
}

template <int dim, typename VectorType>
void
PostProcessorSmoothing<dim, VectorType>::generate_mass_matrix()
{
  dof_handler.distribute_dofs(fe_q);
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  system_matrix->reinit(locally_owned_dofs,
                        locally_owned_dofs,
                        dsp,
                        mpi_communicator);

  QGauss<dim>         quadrature_formula(number_quadrature_points);
  const MappingQ<dim> mapping(
    1, simulation_parameters.fem_parameters.qmapping_all);

  FEValues<dim> fe_values(mapping,
                          fe_q,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const unsigned int dofs_per_cell = fe_q.dofs_per_cell;
<<<<<<< HEAD
  const unsigned int n_q_points    = quadrature_formula.size();

=======
  const unsigned int n_q_points    = number_quadrature_points;

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  AffineConstraints<double> constraints;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

>>>>>>> d040b655 (Created Vorticty and Qcriterion classes)
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  (*system_matrix) = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          local_matrix = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k] = fe_values.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_vf[j] * phi_vf[i]) * fe_values.JxW(q);
                    }
<<<<<<< HEAD
=======
                  // local_rhs(i) += phi_vf[i] * cell_void_fraction *
                  // fe_values_void_fraction.JxW(q);
>>>>>>> d040b655 (Created Vorticty and Qcriterion classes)
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(local_matrix,
                                                 local_dof_indices,
                                                 *system_matrix);
        }
    }
  system_matrix->compress(VectorOperation::add);
<<<<<<< HEAD
}

template <int dim, typename VectorType>
VectorType
PostProcessorSmoothing<dim, VectorType>::solve_L2_projection()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-15;

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs, this->mpi_communicator);

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

  ilu_preconditioner->initialize(*system_matrix, preconditionerOptions);


  solver.solve(*system_matrix,
               completely_distributed_solution,
               system_rhs,
               *ilu_preconditioner);

  constraints.distribute(completely_distributed_solution);

  return completely_distributed_solution;
}

template class PostProcessorSmoothing<2, TrilinosWrappers::MPI::Vector>;
template class PostProcessorSmoothing<3, TrilinosWrappers::MPI::Vector>;
template class PostProcessorSmoothing<2, TrilinosWrappers::MPI::BlockVector>;
template class PostProcessorSmoothing<3, TrilinosWrappers::MPI::BlockVector>;

template <int dim, typename VectorType>
VorticitySmoothing<dim, VectorType>::VorticitySmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points,
  const MPI_Comm &          mpi_communicator)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points,
                                            mpi_communicator)
{}

template <int dim, typename VectorType>
void
VorticitySmoothing<dim, VectorType>::generate_rhs(const VectorType &)
{}

template <int dim, typename VectorType>
void
VorticitySmoothing<dim, VectorType>::calculate_smoothed_field()
=======
  // system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
VorticitySmoothing<dim>::VorticitySmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : PostProcessorSmoothing<dim>(triangulation,
                                simulation_parameters,
                                number_quadrature_points)
{}

template <int dim>
void
VorticitySmoothing<dim>::evaluate_smoothed_field()
{}

template <int dim>
QcriterionSmoothing<dim>::QcriterionSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : PostProcessorSmoothing<dim>(triangulation,
                                simulation_parameters,
                                number_quadrature_points)
{}

template <int dim>
void
QcriterionSmoothing<dim>::evaluate_smoothed_field()
>>>>>>> d040b655 (Created Vorticty and Qcriterion classes)
{}

template class VorticitySmoothing<2, TrilinosWrappers::MPI::Vector>;
template class VorticitySmoothing<3, TrilinosWrappers::MPI::Vector>;
template class VorticitySmoothing<2, TrilinosWrappers::MPI::BlockVector>;
template class VorticitySmoothing<3, TrilinosWrappers::MPI::BlockVector>;

template <int dim, typename VectorType>
QcriterionSmoothing<dim, VectorType>::QcriterionSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points,
  const MPI_Comm &          mpi_communicator)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points,
                                            mpi_communicator)
{}

template <int dim, typename VectorType>
void
QcriterionSmoothing<dim, VectorType>::generate_rhs(const VectorType &solution)
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(
    1, this->simulation_parameters.fem_parameters.qmapping_all);

  const FESystem<dim, dim> fe = this->dof_handler.get_fe();
  FEValues<dim>            fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_q.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  Vector<double>                       local_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);

  this->system_rhs = 0;

  double                           vorticity_on_q_point;
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<2, dim>>      present_velocity_gradients(n_q_points);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          // fe_values[velocities].get_function_gradients(
          // solution, present_velocity_gradients);

          local_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  phi_vf[i] = fe_values.shape_value(i, q);

                  double         p1 = 0.0, r1 = 0.0;
                  Tensor<2, dim> vorticity_vector;
                  Tensor<2, dim> strain_rate_tensor;
                  for (unsigned int j = 0; j < dim; j++)
                    {
                      for (unsigned int k = 0; k < dim; k++)
                        {
                          vorticity_vector[j][k] =
                            0.5 * ((present_velocity_gradients[q][j][k]) -
                                   present_velocity_gradients[q][k][j]);
                          strain_rate_tensor[j][k] =
                            0.5 * ((present_velocity_gradients[q][j][k]) +
                                   present_velocity_gradients[q][k][j]);
                        }
                    }
                  for (unsigned int m = 0; m < dim; m++)
                    {
                      for (unsigned int n = 0; n < dim; n++)
                        {
                          p1 +=
                            (vorticity_vector[m][n]) * (vorticity_vector[m][n]);
                          r1 += (strain_rate_tensor[m][n]) *
                                (strain_rate_tensor[m][n]);
                        }
                    }
                  vorticity_on_q_point = 0.5 * (p1 - r1);

                  local_rhs(i) +=
                    phi_vf[i] * vorticity_on_q_point * fe_values.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          this->constraints.distribute_local_to_global(local_rhs,
                                                       local_dof_indices,
                                                       this->system_rhs);
        }
    }
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim, typename VectorType>
void
QcriterionSmoothing<dim, VectorType>::calculate_smoothed_field()
{}

template class QcriterionSmoothing<2, TrilinosWrappers::MPI::Vector>;
template class QcriterionSmoothing<3, TrilinosWrappers::MPI::Vector>;
template class QcriterionSmoothing<2, TrilinosWrappers::MPI::BlockVector>;
template class QcriterionSmoothing<3, TrilinosWrappers::MPI::BlockVector>;
