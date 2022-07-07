#include <solvers/post_processors_smoothing.h>

template <int dim>
PostProcessorSmoothing<dim>::PostProcessorSmoothing(
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  SimulationParameters<dim> simulation_parameters,
  unsigned int              number_quadrature_points)
  : fe_q(1)
  , dof_handler(*triangulation)
  , simulation_parameters(simulation_parameters)
  , number_quadrature_points(number_quadrature_points)
{}

template <int dim>
void
PostProcessorSmoothing<dim>::generate_mass_matrix()
{
  QGauss<dim> quadrature_formula(number_quadrature_points);
  dof_handler.distribute_dofs(fe_q);
  const MappingQ<dim> mapping(
    1, simulation_parameters.fem_parameters.qmapping_all);

  FEValues<dim> fe_values(mapping,
                          fe_q,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);

  const unsigned int dofs_per_cell = fe_q.dofs_per_cell;
  const unsigned int n_q_points    = number_quadrature_points;

  IndexSet locally_owned_dofs =
    dof_handler.locally_owned_dofs();
  AffineConstraints<double> constraints;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler,
                                          locally_relevant_dofs);
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler,
                                          constraints);
  constraints.close();
  
  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     local_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  phi_vf(dofs_per_cell);

  system_rhs    = 0;
  system_matrix = 0;

  for (const auto &cell :
       dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          local_matrix = 0;
          local_rhs    = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_vf[k]      = fe_values.shape_value(k, q);
                }
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Assemble L2 projection
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_vf[j] * phi_vf[i]) *
                          fe_values.JxW(q);
                    }
                  //local_rhs(i) += phi_vf[i] * cell_void_fraction * fe_values_void_fraction.JxW(q);
                }
            }
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(
            local_matrix,
            local_dof_indices,
            system_matrix);
        }
    }
  system_matrix.compress(VectorOperation::add);
  //system_rhs_void_fraction.compress(VectorOperation::add);
}

template <int dim>
void
PostProcessorSmoothing<dim>::evaluate_smoothed_field()
{}


template class PostProcessorSmoothing<2>;
template class PostProcessorSmoothing<3>;
