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

  // Vector<double>     local_rhs(dofs_per_cell);
}

template <int dim>
void
PostProcessorSmoothing<dim>::evaluate_smoothed_field()
{}


template class PostProcessorSmoothing<2>;
template class PostProcessorSmoothing<3>;
