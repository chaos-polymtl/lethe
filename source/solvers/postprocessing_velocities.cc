#include <solvers/postprocessing_velocities.h>

template TrilinosWrappers::MPI::Vector
AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>::
calculate_average_velocities(
  const TrilinosWrappers::MPI::Vector &local_evaluation_point,
  const Parameters::SimulationControl &simulation_control,
  const Parameters::PostProcessing &   post_processing,
  const IndexSet &                     locally_owned_dofs,
  const MPI_Comm &                     mpi_communicator);

template TrilinosWrappers::MPI::Vector
AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>::
calculate_average_velocities(
  const TrilinosWrappers::MPI::Vector &local_evaluation_point,
  const Parameters::SimulationControl &simulation_control,
  const Parameters::PostProcessing &   post_processing,
  const IndexSet &                     locally_owned_dofs,
  const MPI_Comm &                     mpi_communicator);

/*
template TrilinosWrappers::MPI::BlockVector
AverageVelocities<2, TrilinosWrappers::MPI::BlockVector, std::vector<IndexSet>>::
calculate_average_velocities(
  const TrilinosWrappers::MPI::BlockVector &local_evaluation_point,
  const Parameters::SimulationControl      &simulation_control,
  const Parameters::PostProcessing &        post_processing,
  const std::vector<IndexSet> &             locally_owned_dofs,
  const MPI_Comm &                          mpi_communicator);

template TrilinosWrappers::MPI::BlockVector
AverageVelocities<3, TrilinosWrappers::MPI::BlockVector, std::vector<IndexSet>>::
calculate_average_velocities(
  const TrilinosWrappers::MPI::BlockVector &local_evaluation_point,
  const Parameters::SimulationControl      &simulation_control,
  const Parameters::PostProcessing &        post_processing,
  const std::vector<IndexSet> &             locally_owned_dofs,
  const MPI_Comm &                          mpi_communicator); */