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

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <core/grids.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <dem/dem.h>
#include <dem/dem_properties.h>

#include "core/bdf.h"
#include "core/grids.h"
#include "core/manifolds.h"
#include "core/time_integration_utilities.h"
#include "solvers/gls_navier_stokes.h"



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
  GLSVANSSolver(SimulationParameters<dim> &nsparam);

  ~GLSVANSSolver();

  virtual void
  solve() override;

private:
  void
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  void
  update_solution_and_constraints();

  void
  initialize_void_fraction();

  void
  read_dem();

  void
  calculate_void_fraction(const double time);

  void
  assemble_L2_projection_void_fraction();

  void
  solve_L2_system_void_fraction();

  void
  post_processing();

  virtual void
  iterate() override;

  virtual void
  first_iteration() override;

  /**
   * @brief asocciate the degrees of freedom to each vertex of the finite elements
   * and initialize the void fraction
   */
  virtual void
  setup_dofs() override;


  /**
   * @brief finish_time_step
   * Finishes the time step
   * Post-processing and time stepping
   */
  virtual void
  finish_time_step_fd();

protected:
  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLS();

  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  /**
   * @brief a function for adding data vectors to the data_out object for
   * post_processing additional results
   */
  virtual void
  output_field_hook(DataOut<dim> &data_out) override;

  /**
   *Member Variables
   */

protected:
private:
  DoFHandler<dim> void_fraction_dof_handler;
  FE_Q<dim>       fe_void_fraction;


  MappingQGeneric<dim>                 particle_mapping;
  Particles::ParticleHandler<dim, dim> particle_handler;

  IndexSet locally_owned_dofs_voidfraction;
  IndexSet locally_relevant_dofs_voidfraction;

  // Solution of the void fraction at previous time steps
  TrilinosWrappers::MPI::Vector void_fraction_m1;
  TrilinosWrappers::MPI::Vector void_fraction_m2;
  TrilinosWrappers::MPI::Vector void_fraction_m3;

  TrilinosWrappers::MPI::Vector nodal_void_fraction_relevant;
  TrilinosWrappers::MPI::Vector nodal_void_fraction_owned;

  TrilinosWrappers::SparseMatrix system_matrix_void_fraction;
  TrilinosWrappers::MPI::Vector  system_rhs_void_fraction;
  TrilinosWrappers::SparseMatrix complete_system_matrix_void_fraction;
  TrilinosWrappers::MPI::Vector  complete_system_rhs_void_fraction;
  TrilinosWrappers::SparseMatrix mass_matrix;
  TrilinosWrappers::MPI::Vector  diagonal_of_mass_matrix;
  IndexSet                       active_set;

  // Vectors for drag calculation

  TrilinosWrappers::SparseMatrix system_matrix_drag;
  TrilinosWrappers::MPI::Vector  system_rhs_drag;

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;
  AffineConstraints<double>                          void_fraction_constraints;

  Parameters::SimulationControl::TimeSteppingMethod scheme;

  const bool   PSPG        = true;
  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
};

#endif
