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

#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/time_integration_utilities.h>

#include <solvers/gls_navier_stokes.h>
#include <solvers/postprocessing_cfd.h>

#include <dem/dem.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/vans_assemblers.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>



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
  GLSVANSSolver(CFDDEMSimulationParameters<dim> &nsparam);

  ~GLSVANSSolver();

  virtual void
  solve() override;

private:
  void
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  void
  update_solution_and_constraints();

  void
  particle_centered_method();

  void
  quadrature_centered_sphere_method(bool load_balance_step);

  void
  satellite_point_method();

  void
  solve_L2_system_void_fraction();

  void
  read_dem();

protected:
  /**
   * @brief asocciate the degrees of freedom to each vertex of the finite elements
   * and initialize the void fraction
   */
  virtual void
  setup_dofs() override;

  virtual void
  iterate() override;

  void
  initialize_void_fraction();

  void
  calculate_void_fraction(const double time, bool load_balance_step);

  void
  vertices_cell_mapping();

  virtual void
  monitor_mass_conservation();

  /**
   * @brief finish_time_step
   * Finishes the time step
   * Post-processing and time stepping
   */

  virtual void
  finish_time_step_fd();

  /**
   *  @brief Assembles the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data) override;

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data) override;

  /**
   * @brief sets up the vector of assembler functions
   */
  void
  setup_assemblers() override;


  /**
   * @brief Copy local cell information to global matrix
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copy local cell rhs information to global rhs
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief a function for adding data vectors to the data_out object for
   * post_processing additional results
   */
  virtual void
  output_field_hook(DataOut<dim> &data_out) override;

  void
  percolate_void_fraction();

  /**
   *Member Variables
   */

  CFDDEMSimulationParameters<dim> cfd_dem_simulation_parameters;

  DoFHandler<dim> void_fraction_dof_handler;
  FE_Q<dim>       fe_void_fraction;

  MappingQGeneric<dim> particle_mapping;

  IndexSet locally_owned_dofs_voidfraction;
  IndexSet locally_relevant_dofs_voidfraction;

  // Solution of the void fraction at previous time steps
  std::vector<TrilinosWrappers::MPI::Vector> previous_void_fraction;

  TrilinosWrappers::MPI::Vector nodal_void_fraction_relevant;
  TrilinosWrappers::MPI::Vector nodal_void_fraction_owned;

  // Assemblers for the particle_fluid interactions
  std::vector<std::shared_ptr<ParticleFluidAssemblerBase<dim>>>
    particle_fluid_assemblers;

  TrilinosWrappers::SparseMatrix system_matrix_void_fraction;
  TrilinosWrappers::MPI::Vector  system_rhs_void_fraction;
  TrilinosWrappers::SparseMatrix complete_system_matrix_void_fraction;
  TrilinosWrappers::MPI::Vector  complete_system_rhs_void_fraction;
  TrilinosWrappers::SparseMatrix mass_matrix;
  TrilinosWrappers::MPI::Vector  diagonal_of_mass_matrix;
  IndexSet                       active_set;

  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;
  AffineConstraints<double>                          void_fraction_constraints;

  const bool   PSPG        = true;
  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
  double       pressure_drop;

  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;

protected:
  Particles::ParticleHandler<dim, dim> particle_handler;
};

#endif
