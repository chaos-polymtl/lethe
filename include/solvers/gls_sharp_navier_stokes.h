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
  GLSSharpNavierStokesSolver(SimulationParameters<dim> &nsparam);
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
  // which cell a point falls in vertices_to_cell is a vector of vector of dof
  // handler active cell iterator each element i of the vector is a vector of
  // all the cell in contact with the vertex i.
  void
  vertices_cell_mapping();

  // Defines the particle structure and value based on the parameter file.
  void
  define_particles();

  // Evaluate the forces applied on each of the IB particles.
  void
  force_on_ib();

  // Modify the system matrix to impose IB condition.
  void
  sharp_edge();

  // Write in  the ouput file the force , velocity , position of each of the particles.
  void
  write_force_ib();

  // Integrate the particle velocity and position based on the forces and torque applied to it.
  void
  integrate_particles();

  // Store the solution for the particles dynamics parameters for integration. Defines the table to store the history of each of the particles.
  void
  finish_time_step_particles();

  // Evaluate the L2 error on the computational domain if an analytical solution is given.
  double
  calculate_L2_error_particles();

  virtual void
  postprocess_fd(bool firstIter) override;

  // Allow a refinement around each of the particles.
  void
  refine_ib();

  // Modified version of assemble_matrix_and_rhs to include the presence of extra steps.
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
