/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <solvers/tracer_drift_velocity.h>

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class Parameters::TracerDriftVelocity<2>;
template class Parameters::TracerDriftVelocity<3>;
