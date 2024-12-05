// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_linear_subequations_solver_h
#define lethe_vof_linear_subequations_solver_h

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>
#include <solvers/vof_subequations_interface.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

/**
 * @brief Generalized solver for VOF linear subequations.
 *
 * @tparam dim Number of dimensions of the problem.
 *
 * @ingroup solvers
 */
template <int dim>
class VOFLinearSubequationsSolver : public PhysicsLinearSubequationsSolver
{
public:
  /**
   * @brief Constructor for linear subequations solvers within the VOF
   * auxiliary physics
   *
   * @param[in] p_subequation_id Identifier corresponding to the subequation.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_subequation_verbosity Parameter indicating the verbosity level
   * of the solver.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_multiphysics_interface Multiphysics interface object used to
   * get information from physics.
   *
   * @param[in,out] p_subequations_interface Subequations interface object used
   * to get information from other subequations and store information from the
   * current one.
   */
  VOFLinearSubequationsSolver(
    const VOFSubequationsID         &p_subequation_id,
    const SimulationParameters<dim> &p_simulation_parameters,
    const Parameters::Verbosity     &p_subequation_verbosity,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                  &p_triangulation,
    MultiphysicsInterface<dim>    *p_multiphysics_interface,
    VOFSubequationsInterface<dim> *p_subequations_interface)
    : PhysicsLinearSubequationsSolver(p_pcout)
    , subequation_id(p_subequation_id)
    , subequations_interface(p_subequations_interface)
    , multiphysics_interface(p_multiphysics_interface)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , dof_handler(*triangulation)
    , linear_solver_verbosity(
        p_simulation_parameters.linear_solver.at(PhysicsID::VOF).verbosity)
    , subequation_verbosity(p_subequation_verbosity)
  {}

  /**
   * @brief Default destructor.
   */
  ~VOFLinearSubequationsSolver() = default;

  /**
   * @brief Set up the DofHandler and the degree of freedom associated with
   * the subequation.
   */
  void
  setup_dofs() override;

  /**
   * @brief Assemble and solve linear system when the equation to solve is
   * linear without using the non-linear solver interface.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  void
  solve(const bool &is_post_mesh_adaptation = false) override;

protected:
  /**
   * @brief Assemble the matrix and the right-hand side (rhs) associated with
   * the solver.
   */
  virtual void
  assemble_system_matrix_and_rhs() = 0;

  /**
   * @brief Solve linear system of equation using a strategy appropriate
   * for the partial differential equation.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  void
  solve_linear_system_and_update_solution(
    const bool &is_post_mesh_adaptation = false) override;


  const VOFSubequationsID        subequation_id;
  VOFSubequationsInterface<dim> *subequations_interface;
  MultiphysicsInterface<dim>
    *multiphysics_interface; // to get VOF DoFHandler and solution

  // Parameters
  const SimulationParameters<dim> &simulation_parameters;

  // Core elements
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  DoFHandler<dim>                                              dof_handler;
  std::shared_ptr<FiniteElement<dim>>                          fe;

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>    mapping;
  std::shared_ptr<Quadrature<dim>> cell_quadrature;

  // Solution storage
  IndexSet                                           locally_owned_dofs;
  IndexSet                                           locally_relevant_dofs;
  GlobalVectorType                                   evaluation_point;
  GlobalVectorType                                   present_solution;
  GlobalVectorType                                   system_rhs;
  AffineConstraints<double>                          constraints;
  TrilinosWrappers::SparseMatrix                     system_matrix;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // Verbosity
  const Parameters::Verbosity linear_solver_verbosity;
  const Parameters::Verbosity subequation_verbosity;
};

#endif
