// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_h
#define lethe_fluid_dynamics_vans_h

#include <core/manifolds.h>
#include <core/parameters.h>

#include <solvers/fluid_dynamics_matrix_based.h>
#include <solvers/postprocessing_cfd.h>

#include <dem/dem.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/particle_projector.h>
#include <fem-dem/vans_assemblers.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

using namespace dealii;

/**
 * A solver class for the VANS equation using GLS stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 */

template <int dim>
class FluidDynamicsVANS : public FluidDynamicsMatrixBased<dim>
{
public:
  FluidDynamicsVANS(CFDDEMSimulationParameters<dim> &nsparam);

  ~FluidDynamicsVANS();

  virtual void
  solve() override;

private:
  void
  assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix);

  void
  update_solution_and_constraints();

  void
  read_dem();

protected:
  /**
   * @brief associates the degrees of freedom to each vertex of the finite elements
   * and initialize the void fraction
   */
  virtual void
  setup_dofs() override;

  virtual void
  iterate() override;

  void
  calculate_void_fraction(const double time);

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
   * @brief Assembles the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Assembles the local matrix for a given cell.
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
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief Assembles the local rhs for a given cell
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
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data) override;

  /**
   * @brief sets up the vector of assembler functions
   */
  void
  setup_assemblers() override;


  /**
   * @brief Copies local cell information to global matrix
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief Copies local cell rhs information to global rhs
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data) override;

  /**
   * @brief a function for adding data vectors to the data_out object for
   * post_processing additional results
   */
  virtual std::vector<OutputStruct<dim, GlobalVectorType>>
  get_output_struct_hook() override;


  /**
   * Member Variables
   */

  CFDDEMSimulationParameters<dim> cfd_dem_simulation_parameters;

  MappingQGeneric<dim> particle_mapping;

  // Assemblers for the particle_fluid interactions
  std::vector<std::shared_ptr<ParticleFluidAssemblerBase<dim>>>
    particle_fluid_assemblers;


  const bool   PSPG        = true;
  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
  double       pressure_drop;

  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;

  Particles::ParticleHandler<dim, dim> particle_handler;

  ParticleProjector<dim> particle_projector;

  bool           has_periodic_boundaries;
  Tensor<1, dim> periodic_offset;
  unsigned int   periodic_direction;
};

#endif
