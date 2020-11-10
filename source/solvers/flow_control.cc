#include <solvers/flow_control.h>


template <int dim, typename VectorType>
void
FlowControl<dim, VectorType>::calculate_flow_rate(
  const DoFHandler<dim> &               dof_handler,
  const VectorType &                    present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::FEM &               fem_parameters,
  const MPI_Comm &                      mpi_communicator)
{
  unsigned int boundary_id = flow_control.id_flow_control;

  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const MappingQ<dim>       mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim - 1>           face_quadrature_formula(fe.degree + 1);
  const unsigned int        n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      velocity_values(n_q_points);
  Tensor<1, dim>                   normal_vector;

  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_JxW_values | update_normal_vectors);

  // Resetting next flow rate and area prior new calculation
  flow_rate_n = 0;
  area        = 0;

  // Calculating area and volumetric flow rate at the inlet flow
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              if (cell->face(face)->at_boundary())
                {
                  fe_face_values.reinit(cell, face);
                  if (cell->face(face)->boundary_id() == boundary_id)
                    {
                      for (unsigned int q = 0; q < n_q_points; q++)
                        {
                          area += fe_face_values.JxW(q);
                          normal_vector = fe_face_values.normal_vector(q);
                          fe_face_values[velocities].get_function_values(
                            present_solution, velocity_values);
                          flow_rate_n += velocity_values[q] * normal_vector *
                                         fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }
    }

  area        = Utilities::MPI::sum(area, mpi_communicator);
  flow_rate_n = Utilities::MPI::sum(flow_rate_n, mpi_communicator);
}

template <int dim, typename VectorType>
void
FlowControl<dim, VectorType>::calculate_beta(
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::SimulationControl & simulation_control,
  const unsigned int &                  step_number)
{
  const double dt          = simulation_control.dt;
  const double flow_rate_0 = flow_control.flow_rate;

  beta_n = beta_n1;

  // If flow is reached with no force, adjusted is set to true,
  // threshold factor may be decreased
  if (abs(beta_n - 0) < 1e-6 && step_number > 1)
    {
      if (abs(flow_rate_n) < abs(flow_rate_0))
        adjusted = true;
      else if (abs(flow_rate_1n - flow_rate_n) <
                 abs(flow_rate_n - flow_rate_0) &&
               threshold_factor > 1.00125)
        threshold_factor = 1 + 0.5 * (threshold_factor - 1);
    }

  // If flow rate is over the desired value at time step 1 and beta applied at
  // time step 2 decreased it under the value.
  if (step_number == 3 && abs(flow_rate_n) < abs(flow_rate_0) &&
      abs(flow_rate_1n) > abs(flow_rate_0))
    adjusted = true;

  if (step_number <= 1)
    {
      // Initial beta
      beta_n1 = flow_control.beta_0;
    }
  else if (abs(flow_rate_n) > abs(flow_rate_0) &&
           abs(flow_rate_n) < threshold_factor * abs(flow_rate_0) &&
           adjusted == false)
    {
      // If the flow rate is between targeted flow rate value and the threshold
      // and if it didn't reached it the value, it decreases by itself
      // (no force applied)
      beta_n1 = 0;
    }
  else if (step_number == 2)
    {
      // 2nd time step if not in the
      beta_n1 = 0.5 * (flow_rate_n - flow_rate_0) / (area * dt);
    }
  else
    {
      // Standard flow controller.
      // Calculate the new beta to control the flow.
      // If desired flow rate is reached, new beta only maintains the force to
      // keep the flow at the desired value. Is so, if calculated beta is
      // negative it's set to 0 to avoided +/- oscillations
      beta_n1 =
        beta_n - (flow_rate_0 - 2 * flow_rate_n + flow_rate_1n) / (area * dt);
      // If flow rate was already adjusted and beta is the same sign
      // than flow_rate_0
      if (flow_rate_0 * beta_n1 > 0 && adjusted == true)
        beta_n1 = 0;
    }

  flow_rate_1n = flow_rate_n;

  // Setting beta coefficient only to new time step
  if (flow_control.flow_direction == 0)
    beta[0] = beta_n1; // beta = f_x
  else if (flow_control.flow_direction == 1)
    beta[1] = beta_n1; // beta = f_y
  else if (flow_control.flow_direction == 2)
    beta[2] = beta_n1; // beta = f_z
}


template <int dim, typename VectorType>
double
FlowControl<dim, VectorType>::get_flow_rate(
  const DoFHandler<dim> &               dof_handler,
  const VectorType &                    present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::FEM &               fem_parameters,
  const MPI_Comm &                      mpi_communicator)
{
  calculate_flow_rate(dof_handler,
                      present_solution,
                      flow_control,
                      fem_parameters,
                      mpi_communicator);
  return flow_rate_n;
}


template <int dim, typename VectorType>
Tensor<1, dim>
FlowControl<dim, VectorType>::get_beta(
  const DoFHandler<dim> &               dof_handler,
  const VectorType &                    present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::SimulationControl & simulation_control,
  const Parameters::FEM &               fem_parameters,
  const unsigned int &                  step_number,
  const MPI_Comm &                      mpi_communicator)
{
  calculate_flow_rate(dof_handler,
                      present_solution,
                      flow_control,
                      fem_parameters,
                      mpi_communicator);

  calculate_beta(flow_control, simulation_control, step_number);

  return beta;
}


template class FlowControl<2, TrilinosWrappers::MPI::Vector>;
template class FlowControl<3, TrilinosWrappers::MPI::Vector>;
template class FlowControl<2, TrilinosWrappers::MPI::BlockVector>;
template class FlowControl<3, TrilinosWrappers::MPI::BlockVector>;
