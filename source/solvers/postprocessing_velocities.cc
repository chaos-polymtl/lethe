#include <solvers/postprocessing_velocities.h>


template std::vector<Tensor<1, 2>>
PostprocessingVelocities<2, TrilinosWrappers::MPI::Vector>::
  calculate_velocity_fluctuations(
    const DoFHandler<2> &                dof_handler,
    const TrilinosWrappers::MPI::Vector &evaluation_point,
    const Parameters::SimulationControl &simulation_control,
    const Parameters::FEM &              fem_parameters,
    const double &                       step_number,
    const double &                       bulk_velocity,
    const MPI_Comm &                     mpi_communicator);

template std::vector<Tensor<1, 3>>
PostprocessingVelocities<3, TrilinosWrappers::MPI::Vector>::
  calculate_velocity_fluctuations(
    const DoFHandler<3> &                dof_handler,
    const TrilinosWrappers::MPI::Vector &evaluation_point,
    const Parameters::SimulationControl &simulation_control,
    const Parameters::FEM &              fem_parameters,
    const double &                       step_number,
    const double &                       bulk_velocity,
    const MPI_Comm &                     mpi_communicator);

template std::vector<Tensor<1, 2>>
PostprocessingVelocities<2, TrilinosWrappers::MPI::BlockVector>::
  calculate_velocity_fluctuations(
    const DoFHandler<2> &                     dof_handler,
    const TrilinosWrappers::MPI::BlockVector &evaluation_point,
    const Parameters::SimulationControl &     simulation_control,
    const Parameters::FEM &                   fem_parameters,
    const double &                            step_number,
    const double &                            bulk_velocity,
    const MPI_Comm &                          mpi_communicator);

template std::vector<Tensor<1, 3>>
PostprocessingVelocities<3, TrilinosWrappers::MPI::BlockVector>::
  calculate_velocity_fluctuations(
    const DoFHandler<3> &                     dof_handler,
    const TrilinosWrappers::MPI::BlockVector &evaluation_point,
    const Parameters::SimulationControl &     simulation_control,
    const Parameters::FEM &                   fem_parameters,
    const double &                            step_number,
    const double &                            bulk_velocity,
    const MPI_Comm &                          mpi_communicator);
