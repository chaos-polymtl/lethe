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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef LETHE_GLSSHARPNS_H
#define LETHE_GLSSHARPNS_H

#include <core/ib_particle.h>
#include <solvers/gls_navier_stokes.h>

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization and
 * Sharp-Edge immersed boundaries
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Lucka Barbeau, Bruno Blais, 2020
 */

template <int dim>
class GLSSharpNavierStokesSolver : public GLSNavierStokesSolver<dim>
{
public:
  GLSSharpNavierStokesSolver(NavierStokesSolverParameters<dim> &nsparam,
                             const unsigned int                 degreeVelocity,
                             const unsigned int                 degreePressure);
  ~GLSSharpNavierStokesSolver();

  void
  solve();

private:
  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLS();


  // BB - TODO This explanation needs to be made clearer. Adjacent, Adjacent_2
  // and Adjacent_3 needs to be renamed if possible to a clearer notation

  // Map the vertex index to the cell that include that vertex used later in
  // which cell a point falls in vertices_to_cell is a vector of vectof of dof
  // handler active cell iterator each element i of the vector is a vector of
  // all the cell in contact with the vertex i
  void
  vertices_cell_mapping();

  // BB - TODO The particles structure should be refactored to use a small class
  // to store the information or a struct instead of just using a vector where
  // the things are hardcoded within.
  void
  define_particles();

  void
  force_on_ib();

  void
  sharp_edge();

  void
  write_force_ib();



  double
  calculate_L2_error_particles();

  virtual void
  postprocess(bool firstIter) override;

  void
  refine_ib();

  void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  /**
   * Members
   */
private:
  std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator>>
                               vertices_to_cell;
  const bool                   SUPG        = false;
  const bool                   PSPG        = true;
  const double                 GLS_u_scale = 1;
  std::vector<IBParticle<dim>> particles;


  std::vector<TableHandler> table_f;
  std::vector<TableHandler> table_t;
};


#endif
