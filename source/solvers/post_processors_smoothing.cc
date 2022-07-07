#include <solvers/post_processors_smoothing.h>

template <int dim, typename VectorType>
PostProcessorSmoothing<dim, VectorType>::PostProcessorSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : fe_q(1)
  , dof_handler(*triangulation)
  , simulation_parameters(simulation_parameters)
  , number_quadrature_points(number_quadrature_points)
{}

template <int dim, typename VectorType>
void
PostProcessorSmoothing<dim, VectorType>::generate_mass_matrix()
{
  QGauss<dim> quadrature_formula(number_quadrature_points);
  dof_handler.distribute_dofs(fe_q);
  const MappingQ<dim> mapping(
    1, simulation_parameters.fem_parameters.qmapping_all);

  FEValues<dim> fe_values(mapping,
                          fe_q,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  const unsigned int dofs_per_cell = fe_q.dofs_per_cell;
  const unsigned int n_q_points    = number_quadrature_points;

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  AffineConstraints<double> constraints;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  system_matrix = 0;

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
VorticitySmoothing<dim, VectorType>::VorticitySmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points)
{}

template <int dim, typename VectorType>
void
VorticitySmoothing<dim, VectorType>::evaluate_smoothed_field(const VectorType &)
{}

template <int dim, typename VectorType>
QcriterionSmoothing<dim, VectorType>::QcriterionSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : PostProcessorSmoothing<dim, VectorType>(triangulation,
                                            simulation_parameters,
                                            number_quadrature_points)
{}

template <int dim, typename VectorType>
void
QcriterionSmoothing<dim, VectorType>::evaluate_smoothed_field(
  const VectorType &solution)
{
  QGauss<dim>         quadrature_formula(this->number_quadrature_points);
  const MappingQ<dim> mapping(
    1, this->simulation_parameters.fem_parameters.qmapping_all);

  FEValues<dim> fe_values(mapping,
                          this->fe_q,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = this->fe_q.dofs_per_cell;
  const unsigned int n_q_points    = this->number_quadrature_points;

  IndexSet locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  AffineConstraints<double> constraints;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          locally_relevant_dofs);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->dof_handler, constraints);
  constraints.close();

  Vector<double>                       local_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  system_rhs = 0;

  double                           vorticity_on_q_point;
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<2, dim>>      present_velocity_gradients(n_q_points);

  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[velocities].get_function_gradients(
            solution, present_velocity_gradients);

          local_rhs = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  phi_vf[i] = fe_values.shape_value(i, q);

                  double                      p1 = 0.0, r1 = 0.0;
                  std::vector<Tensor<2, dim>> vorticity_vector(
                    present_velocity_gradients.size());
                  std::vector<Tensor<2, dim>> strain_rate_tensor(
                    present_velocity_gradients.size());
                  for (unsigned int j = 0; j < dim; j++)
                    {
                      for (unsigned int k = 0; k < dim; k++)
                        {
                          vorticity_vector[q][j][k] =
                            0.5 * ((present_velocity_gradients[q][j][k]) -
                                   present_velocity_gradients[q][k][j]);
                          strain_rate_tensor[q][j][k] =
                            0.5 * ((present_velocity_gradients[q][j][k]) +
                                   present_velocity_gradients[q][k][j]);
                        }
                    }
                  for (unsigned int m = 0; m < dim; m++)
                    {
                      for (unsigned int n = 0; n < dim; n++)
                        {
                          p1 += (vorticity_vector[q][m][n]) *
                                (vorticity_vector[q][m][n]);
                          r1 += (strain_rate_tensor[q][m][n]) *
                                (strain_rate_tensor[q][m][n]);
                        }
                    }
                  vorticity_on_q_point = 0.5 * (p1 - r1);

                  local_rhs(i) +=
                    phi_vf[i] * vorticity_on_q_point * fe_values.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(local_rhs,
                                                 local_dof_indices,
                                                 system_rhs);
        }
    }
  system_rhs.compress(VectorOperation::add);
}


template class PostProcessorSmoothing<2, TrilinosWrappers::MPI::Vector>;
template class PostProcessorSmoothing<3, TrilinosWrappers::MPI::Vector>;
template class PostProcessorSmoothing<2, TrilinosWrappers::MPI::BlockVector>;
template class PostProcessorSmoothing<3, TrilinosWrappers::MPI::BlockVector>;