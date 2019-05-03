#include "glsNS.h"

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is valid before we actually compile the solver
// This greatly helps with debugging
extern template class GLSNavierStokesSolver<2>;
extern template class GLSNavierStokesSolver<3>;
