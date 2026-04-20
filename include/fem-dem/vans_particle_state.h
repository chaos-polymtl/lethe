// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vans_particle_state_h
#define lethe_vans_particle_state_h

#include <core/boundary_conditions.h>
#include <core/parameters.h>

#include <fem-dem/particle_projector.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Bundle of particle-side state shared by the VANS fluid solvers.
 *
 * Holds the particle mapping, particle handler, particle projector, and the
 * periodic-boundary bookkeeping that both FluidDynamicsVANS and
 * FluidDynamicsVANSMatrixFree need. The members are public by design: the
 * VANS solvers use them directly, exactly as if they were their own fields.
 *
 * @tparam dim Spatial dimension of the simulation (2 or 3)
 */
template <int dim>
class VANSParticleState
{
public:
  /**
   * @brief Construct the particle-side state used by a VANS solver.
   *
   * @param[in] triangulation The triangulation used for the CFD simulation.
   * @param[in] void_fraction_parameters Parameters controlling the void
   * fraction computation.
   * @param[in] linear_solver_parameters Linear solver parameters used by the
   * particle projector for the L2 projection of the void fraction.
   * @param[in] void_fraction_order FE degree used to interpolate the void
   * fraction.
   * @param[in] simplex True if the mesh uses simplex elements.
   * @param[in] pcout Conditional output stream used for logging.
   */
  VANSParticleState(parallel::DistributedTriangulationBase<dim> *triangulation,
                    std::shared_ptr<Parameters::VoidFractionParameters<dim>>
                                                    void_fraction_parameters,
                    const Parameters::LinearSolver &linear_solver_parameters,
                    const unsigned int              void_fraction_order,
                    const bool                      simplex,
                    const ConditionalOStream       &pcout);

  /**
   * @brief Scan boundary conditions and update @p has_periodic_boundaries.
   *
   * Sets @p has_periodic_boundaries to true when exactly one periodic
   * boundary pair is declared. Throws if more than one periodic boundary
   * condition is present, matching the historical VANS restriction.
   *
   * @param[in] boundary_conditions Navier-Stokes boundary conditions.
   */
  void
  scan_periodic_boundaries(
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions);

  /**
   * @brief Set up DoFs and constraints of the particle projector.
   *
   * Calls ParticleProjector::setup_dofs followed by setup_constraints.
   *
   * @param[in] boundary_conditions Navier-Stokes boundary conditions passed
   * through to the projector's constraint setup.
   */
  void
  setup_dofs(
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions);

  /**
   * @brief Load a DEM checkpoint and populate @p particle_handler.
   *
   * Reads the DEM particle handler and triangulation from the files identified
   * by @p dem_file_name and converts the DEM properties into CFD-DEM
   * properties.
   *
   * @param[in] dem_file_name Checkpoint file prefix written by the DEM solver.
   * @param[in] require_3d When true, assert that dim == 3 before loading.
   * Used by the matrix-free VANS solver which only supports 3D.
   * @param[in] exchange_ghosts When true, call exchange_ghost_particles on
   * the handler after loading. Used by the matrix-free VANS solver.
   */
  void
  read_dem(const std::string &dem_file_name,
           const bool         require_3d      = false,
           const bool         exchange_ghosts = false);

  /// Non-owning pointer to the CFD triangulation. Needed to access the
  /// distributed triangulation as non-const for checkpoint loading, since
  /// ParticleHandler::get_triangulation returns a const reference.
  parallel::DistributedTriangulationBase<dim> *triangulation;

  /// Mapping used by the particle handler.
  MappingQGeneric<dim> particle_mapping;

  /// Particle handler owning the CFD-DEM particles.
  Particles::ParticleHandler<dim, dim> particle_handler;

  /// Particle projector computing the void fraction and particle fields.
  ParticleProjector<dim> particle_projector;

  /// True if the simulation has periodic boundary conditions.
  bool has_periodic_boundaries;

  /// Offset vector between the two periodic boundaries (unused today; kept
  /// for parity with the prior duplicated state).
  Tensor<1, dim> periodic_offset;

  /// Direction index of the periodic boundary pair (unused today; kept for
  /// parity with the prior duplicated state).
  unsigned int periodic_direction;
};

#endif
