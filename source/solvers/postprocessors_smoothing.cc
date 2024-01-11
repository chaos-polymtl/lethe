#include <solvers/postprocessors_smoothing.h>

template <int dim, typename VectorType>
PostProcessorSmoothing<dim, VectorType>::PostProcessorSmoothing(
  const parallel::DistributedTriangulationBase<dim> &triangulation,
  const SimulationParameters<dim>                   &simulation_parameters,
  const unsigned int                                &number_quadrature_points)
  : fe_q(1)
  , dof_handler(triangulation)
  , simulation_parameters(simulation_parameters)
  , number_quadrature_points(number_quadrature_points)
  , mpi_communicator(triangulation.get_communicator())
{
  mapping = std::make_shared<MappingQ<dim>>(
    1, this->simulation_parameters.fem_parameters.qmapping_all);
  system_matrix = std::make_shared<TrilinosWrappers::SparseMatrix>();
}

template <int dim, typename VectorType>
void
PostProcessorSmoothing<dim, VectorType>::generate_mass_matrix()
{
  dof_handler.distribute_dofs(fe_q);
  this->locally_owned_dofs = dof_handler.locally_owned_dofs();

  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             this->locally_owned_dofs,
                                             this->mpi_communicator,
                                             this->locally_relevant_dofs);

  system_matrix->reinit(this->locally_owned_dofs,
                        this->locally_owned_dofs,
                        dsp,
                        this->mpi_communicator);

  QGauss<dim> quadrature_formula(number_quadrature_points);

  FEValues<dim> fe_values(*this->mapping,
                          fe_q,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_q.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_vf(dofs_per_cell);


  (*system_matrix) = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          double h = cell->diameter();

          local_matrix = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values.JxW(q);
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values.shape_value(k, q);
                  grad_phi_vf[k] = fe_values.shape_grad(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection with smoothing
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_vf[j] * phi_vf[i] +
                         h * h * grad_phi_vf[j] * grad_phi_vf[i]) *
                        JxW;
                    }
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(local_matrix,
                                                 local_dof_indices,
                                                 *system_matrix);
        }
    }
  system_matrix->compress(VectorOperation::add);
}

template <int dim, typename VectorType>
const GlobalVectorType &
PostProcessorSmoothing<dim, VectorType>::solve_L2_projection()
{
  // Solve the L2 projection system
  const double linear_solver_tolerance = 1e-12;

  completely_distributed_solution =
    GlobalVectorType(this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverCG solver(solver_control);

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const double ilu_fill = 0.;
  const double ilu_atol = 1e-15;
  const double ilu_rtol = 1.;

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

template <int dim, typename VectorType>
const GlobalVectorType &
PostProcessorSmoothing<dim, VectorType>::calculate_smoothed_field(
  const VectorType             &solution,
  const DoFHandler<dim>        &dof_handler_velocity,
  std::shared_ptr<Mapping<dim>> mapping_velocity)
{
  generate_mass_matrix();
  generate_rhs(solution, dof_handler_velocity, mapping_velocity);
  return solve_L2_projection();
}

template <int dim, typename VectorType>
const DoFHandler<dim> &
PostProcessorSmoothing<dim, VectorType>::get_dof_handler() const
{
  return dof_handler;
}

template class PostProcessorSmoothing<2, GlobalVectorType>;
template class PostProcessorSmoothing<3, GlobalVectorType>;
template class PostProcessorSmoothing<2, GlobalBlockVectorType>;
template class PostProcessorSmoothing<3, GlobalBlockVectorType>;


template <int dim, typename VectorType>
QcriterionPostProcessorSmoothing<dim, VectorType>::
  QcriterionPostProcessorSmoothing(
    const parallel::DistributedTriangulationBase<dim> &triangulation,
    const SimulationParameters<dim>                   &simulation_parameters,
    const unsigned int                                &number_quadrature_points)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points)
{}

template <int dim, typename VectorType>
void
QcriterionPostProcessorSmoothing<dim, VectorType>::generate_rhs(
  const VectorType             &solution,
  const DoFHandler<dim>        &dof_handler_velocity,
  std::shared_ptr<Mapping<dim>> mapping_velocity)
{
  QGauss<dim> quadrature_formula(this->number_quadrature_points);

  const FESystem<dim, dim> fe = this->dof_handler.get_fe();
  FEValues<dim>            fe_values(*this->mapping,
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

  // Fluid dynamics information
  const FESystem<dim, dim> fe_velocity = dof_handler_velocity.get_fe();
  FEValues<dim>            fe_values_velocity(*mapping_velocity,
                                   fe_velocity,
                                   quadrature_formula,
                                   update_quadrature_points | update_gradients);

  for (const auto &cell : dof_handler_velocity.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Because we are looping over the active cells of the
          // dof_handler_velocity.  We need to explicitely build a cell iterator
          // to the same cell, but for the smoother dof_handler
          const auto &smoother_cell =
            typename DoFHandler<dim>::cell_iterator(*cell, &this->dof_handler);
          fe_values.reinit(smoother_cell);

          fe_values_velocity.reinit(cell);

          fe_values_velocity[velocities].get_function_gradients(
            solution, present_velocity_gradients);

          local_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values.JxW(q);
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

                  local_rhs(i) += phi_vf[i] * vorticity_on_q_point * JxW;
                }
            }
          // Associate cell of fluid dof_handler to current (qcriterion)
          // dof_handler
          smoother_cell->get_dof_indices(local_dof_indices);
          this->constraints.distribute_local_to_global(local_rhs,
                                                       local_dof_indices,
                                                       this->system_rhs);
        }
    }
  this->system_rhs.compress(VectorOperation::add);
}

template class QcriterionPostProcessorSmoothing<2, GlobalVectorType>;
template class QcriterionPostProcessorSmoothing<3, GlobalVectorType>;
template class QcriterionPostProcessorSmoothing<2, GlobalBlockVectorType>;
template class QcriterionPostProcessorSmoothing<3, GlobalBlockVectorType>;

template <int dim, typename VectorType>
ContinuityPostProcessorSmoothing<dim, VectorType>::
  ContinuityPostProcessorSmoothing(
    const parallel::DistributedTriangulationBase<dim> &triangulation,
    const SimulationParameters<dim>                   &simulation_parameters,
    const unsigned int                                &number_quadrature_points)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points)
{}

template <int dim, typename VectorType>
void
ContinuityPostProcessorSmoothing<dim, VectorType>::generate_rhs(
  const VectorType             &solution,
  const DoFHandler<dim>        &dof_handler_velocity,
  std::shared_ptr<Mapping<dim>> mapping_velocity)
{
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs = 0;

  // Smoother management
  QGauss<dim> quadrature_formula(this->number_quadrature_points);

  const FESystem<dim, dim> fe = this->dof_handler.get_fe();
  FEValues<dim>            fe_values(*this->mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const unsigned int dofs_per_cell = this->fe_q.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  Vector<double>                       local_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  // Velocity information
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<2, dim>>      present_velocity_gradients(n_q_points);

  const FESystem<dim, dim> fe_velocity = dof_handler_velocity.get_fe();
  FEValues<dim>            fe_values_velocity(*mapping_velocity,
                                   fe_velocity,
                                   quadrature_formula,
                                   update_quadrature_points | update_gradients);

  for (const auto &cell : dof_handler_velocity.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Because we are looping over the active cells of the
          // dof_handler_velocity.  We need to explicitely build a cell iterator
          // to the same cell, but for the smoother dof_handler
          const auto &smoother_cell =
            typename DoFHandler<dim>::cell_iterator(*cell, &this->dof_handler);
          fe_values.reinit(smoother_cell);

          fe_values_velocity.reinit(cell);

          fe_values_velocity[velocities].get_function_gradients(
            solution, present_velocity_gradients);

          local_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values.JxW(q);
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  phi_vf[i]         = fe_values.shape_value(i, q);
                  double continuity = 0;
                  for (unsigned int d = 0; d < dim; d++)
                    {
                      continuity += present_velocity_gradients[q][d][d];
                    }

                  local_rhs(i) += phi_vf[i] * continuity * JxW;
                }
            }
          // Associate cell of velocity dof_handler to current (qcriterion)
          // dof_handler
          smoother_cell->get_dof_indices(local_dof_indices);
          this->constraints.distribute_local_to_global(local_rhs,
                                                       local_dof_indices,
                                                       this->system_rhs);
        }
    }
  this->system_rhs.compress(VectorOperation::add);
}

template class ContinuityPostProcessorSmoothing<2, GlobalVectorType>;
template class ContinuityPostProcessorSmoothing<3, GlobalVectorType>;
template class ContinuityPostProcessorSmoothing<2, GlobalBlockVectorType>;
template class ContinuityPostProcessorSmoothing<3, GlobalBlockVectorType>;
