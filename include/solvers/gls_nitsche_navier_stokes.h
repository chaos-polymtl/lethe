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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020-
 */

#ifndef lethe_gls_nitsche_navier_stokes_h
#define lethe_gls_nitsche_navier_stokes_h

#include "gls_navier_stokes.h"
#include "solid_base.h"

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization
 * and Nitsche immersed boundary method
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim, int spacedim = dim>
class GLSNitscheNavierStokesSolver
  : public GLSNavierStokesSolver<spacedim>
{
public:
  GLSNitscheNavierStokesSolver(NavierStokesSolverParameters<spacedim> &nsparam,
                               const unsigned int                      degreeVelocity,
                               const unsigned int                      degreePressure);
  ~GLSNitscheNavierStokesSolver();

private:

  SolidBase<dim,spacedim>  solid;

  void
  assemble_nitsche_restriction();

  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;
};


#endif
