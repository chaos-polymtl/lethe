#include <solvers/postprocessing_velocities.h>

template class
AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class
AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class
AverageVelocities<2, TrilinosWrappers::MPI::BlockVector,
                    std::vector<IndexSet>>;

template class
AverageVelocities<3, TrilinosWrappers::MPI::BlockVector,
                    std::vector<IndexSet>>;
