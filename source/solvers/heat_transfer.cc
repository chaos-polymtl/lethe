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
{}


template <int dim>
void
HeatTransfer<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{}

template <int dim>
void
HeatTransfer<dim>::finish_simulation()
{}

template <int dim>
void
HeatTransfer<dim>::finish_time_step()
{}

template <int dim>
void
HeatTransfer<dim>::postprocess()
{}

template <int dim>
void
HeatTransfer<dim>::setup_dofs()
{}

template <int dim>
void
HeatTransfer<dim>::set_initial_conditions()
{}

template <int dim>
void
HeatTransfer<dim>::solve_linear_system(const bool initial_step,
                                       const bool renewed_matrix)
{}



template class HeatTransfer<2>;
template class HeatTransfer<3>;
