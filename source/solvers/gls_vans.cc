#include "solvers/gls_vans.h"
// Constructor for class GLS_VANS
template <int dim>
GLSVANSSolver<dim>::GLSVANSSolver(NavierStokesSolverParameters<dim> &p_nsparam,
                                  const unsigned int p_degree_velocity,
                                  const unsigned int p_degree_pressure)
  : GLSNavierStokesSolver<dim>(p_nsparam, p_degree_velocity, p_degree_pressure)

{}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class GLSVANSSolver<2>;
