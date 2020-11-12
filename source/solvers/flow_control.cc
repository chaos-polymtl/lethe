#include <solvers/flow_control.h>

template Tensor<1, 2>
FlowControl<2, TrilinosWrappers::MPI::Vector>::get_beta(
  const DoFHandler<2> &                 dof_handler,
  const TrilinosWrappers::MPI::Vector & present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::SimulationControl & simulation_control,
  const Parameters::FEM &               fem_parameters,
  const double &                        step_number,
  const MPI_Comm &                      mpi_communicator);

template Tensor<1, 3>
FlowControl<3, TrilinosWrappers::MPI::Vector>::get_beta(
  const DoFHandler<3> &                 dof_handler,
  const TrilinosWrappers::MPI::Vector & present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::SimulationControl & simulation_control,
  const Parameters::FEM &               fem_parameters,
  const double &                        step_number,
  const MPI_Comm &                      mpi_communicator);

template Tensor<1, 2>
FlowControl<2, TrilinosWrappers::MPI::BlockVector>::get_beta(
  const DoFHandler<2> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &present_solution,
  const Parameters::DynamicFlowControl &    flow_control,
  const Parameters::SimulationControl &     simulation_control,
  const Parameters::FEM &                   fem_parameters,
  const double &                            step_number,
  const MPI_Comm &                          mpi_communicator);

template Tensor<1, 3>
FlowControl<3, TrilinosWrappers::MPI::BlockVector>::get_beta(
  const DoFHandler<3> &                     dof_handler,
  const TrilinosWrappers::MPI::BlockVector &present_solution,
  const Parameters::DynamicFlowControl &    flow_control,
  const Parameters::SimulationControl &     simulation_control,
  const Parameters::FEM &                   fem_parameters,
  const double &                            step_number,
  const MPI_Comm &                          mpi_communicator);