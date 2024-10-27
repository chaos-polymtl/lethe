// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_phase_gradient_projection_h
#define lethe_vof_phase_gradient_projection_h

#include <solvers/vof_linear_subequations_solver.h>

/**
 * @brief VOF phase fraction gradient L2 projection solver.
 *
 * @tparam dim Number of dimensions of the problem.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 */
template <int dim, typename ScratchDataType>
class VOFPhaseGradientProjection
  : public VOFLinearSubequationsSolver<dim, ScratchDataType>
{
public:
  /**
   * @brief Constructor for the L2 projection of the VOF phase fraction gradient
   * (pfg).
   *
   * @param[in,out] p_subequations Subequations interface object used to get
   * information from other subequations and store information from the current
   * one.
   *
   * @param[in] p_multiphysics Multiphysics interface object used to get
   * information from physics.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_simulation_control Object responsible for the control of
   * steady-state and transient simulations. Contains all the information
   * related to time stepping and the stopping criteria.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   */
  VOFPhaseGradientProjection(
    VOFSubequationsInterface<dim>   *p_subequations,
    MultiphysicsInterface<dim>      *p_multiphysics,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       &p_triangulation,
    std::shared_ptr<SimulationControl> &p_simulation_control,
    const ConditionalOStream           &p_pcout)
    : VOFLinearSubequationsSolver<dim, ScratchDataType>(
        VOFSubequationsID::phase_gradient_projection,
        p_subequations,
        p_multiphysics,
        p_simulation_parameters,
        p_triangulation,
        p_simulation_control,
        p_simulation_parameters.multiphysics.vof_parameters
          .surface_tension_force.verbosity,
        p_pcout)
  {}

  /**
   * @brief Default destructor.
   */
  ~VOFPhaseGradientProjection() = default;

private:
  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * @param[in] cell Cell for which the local matrix is assembled.
   *
   * @param[in] scratch_data Stores the calculated finite element information at
   * Gauss points.
   *
   * @param[in,out] copy_data Stores the results of the assembly over a cell.
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchDataType                                      &scratch_data,
    StabilizedMethodsCopyData                            &copy_data) override;

  /**
   * @brief Assemble the local right-hand side (rhs) for a given cell.
   *
   * @param[in] cell Cell for which the local rhs is assembled.
   *
   * @param[in] scratch_data Stores the calculated finite element information at
   * Gauss points.
   *
   * @param[in,out] copy_data Stores the results of the assembly over a cell.
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchDataType                                      &scratch_data,
    StabilizedMethodsCopyData                            &copy_data) override;
};

#endif
