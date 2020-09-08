/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Toni EL Geitani, Polytechnique Montreal, 2020-
 */

#ifndef lethe_gls_vans_h
#define lethe_gls_vans_h

#include "gls_navier_stokes.h"

using namespace dealii;

/**
 * A solver class for the VANS equation using GLS stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Toni EL Geitani, 2020
 */

template <int dim>
class GLSVANSSolver : public GLSNavierStokesSolver<dim>
{
public:
  GLSVANSSolver(NavierStokesSolverParameters<dim> &nsparam,
                               const unsigned int degree_velocity,
                               const unsigned int degree_pressure);


};

#endif


