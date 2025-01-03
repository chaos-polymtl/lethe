#include <dem/integrator.h>

template class Integrator<2, DEM::SolverType::dem>;
template class Integrator<2, DEM::SolverType::cfd_dem>;
template class Integrator<3, DEM::SolverType::dem>;
template class Integrator<3, DEM::SolverType::cfd_dem>;
