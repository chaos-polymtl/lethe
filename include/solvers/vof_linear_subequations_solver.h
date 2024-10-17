// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_linear_subequations_solver_h
#define lethe_vof_linear_subequations_solver_h

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>

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

#include "vof_subequations_interface.h"

/**
 * @brief Generalized solver for VOF linear subequations.
 *
 * @tparam dim Number of dimensions of the problem.
 *
 * @tparam ScratchDataType Type of scratch data object used for linear system
 * assembly.
 *
 * @ingroup solvers
 */
template <int dim, typename ScratchDataType>
class VOFLinearSubequationsSolver : public PhysicsLinearSubequationsSolver
{
public:
  /**
   * @brief Constructor for linear subequations solvers within the VOF
   * auxiliary physics
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
  VOFLinearSubequationsSolver(
    VOFSubequationsID                p_subequation_id,
    VOFSubequationsInterface<dim>   *p_subequations,
    MultiphysicsInterface<dim>      *p_multiphysics,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       &p_triangulation,
    std::shared_ptr<SimulationControl> &p_simulation_control,
    const Parameters::Verbosity        &p_subequation_verbosity,
    const ConditionalOStream           &p_pcout)
    : PhysicsLinearSubequationsSolver(p_pcout)
    , subequation_id(p_subequation_id)
    , subequations(p_subequations)
    , multiphysics(p_multiphysics)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
    , linear_solver_verbosity(
        p_simulation_parameters.linear_solver.at(PhysicsID::VOF).verbosity)
    , subequation_verbosity(p_subequation_verbosity)
  {
    if (this->simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        const FE_SimplexP<dim> phase_gradient_fe(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->fe      = std::make_shared<FESystem<dim>>(phase_gradient_fe, dim);
        this->mapping = std::make_shared<MappingFE<dim>>(*this->fe);

        this->cell_quadrature =
          std::make_shared<QGaussSimplex<dim>>(this->fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        const FE_Q<dim> phase_gradient_fe(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->fe      = std::make_shared<FESystem<dim>>(phase_gradient_fe, dim);
        this->mapping = std::make_shared<MappingQ<dim>>(this->fe->degree);

        this->cell_quadrature =
          std::make_shared<QGauss<dim>>(this->fe->degree + 1);
      }
  }

  /**
   * @brief Set up the DofHandler and the degree of freedom associated with
   * the physics.
   *
   * @param[in] subequation_id Identifier corresponding to the subequation, for
   * verbosity purpose.
   */
  void
  setup_dofs() override;

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
   * @brief Assemble the matrix associated with the solver
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
   * @param[in] cell Cell for which the local matrix is assembled.
   *
   * @param[in] scratch_data Stores the calculated finite element information at
   * Gauss points.
   *
   * @param[out] copy_data Stores the results of the assembly over a cell.
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchDataType                                      &scratch_data,
    StabilizedMethodsCopyData                            &copy_data) = 0;

  /**
   * @brief Assemble the local right-hand side (rhs) for a given cell.
   *
   * @param[in] cell Cell for which the local rhs is assembled.
   *
   * @param[in] scratch_data Stores the calculated finite element information at
   * Gauss points.
   *
   * @param[out] copy_data Stores the results of the assembly over a cell.
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    ScratchDataType                                      &scratch_data,
    StabilizedMethodsCopyData                            &copy_data) = 0;

  /**
   * @brief Copy local cell matrix information to global matrix.
   *
   * @param[in] copy_data Stores the results of the assembly over a cell.
   */
  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Copy local right-hand side (rhs) information to global rhs.
   *
   * @param[in] copy_data Stores the results of the assembly over a cell.
   */
  void
  copy_local_rhs_to_global_rhs(const StabilizedMethodsCopyData &copy_data);


  VOFSubequationsID              subequation_id;
  VOFSubequationsInterface<dim> *subequations;
  MultiphysicsInterface<dim>
    *multiphysics; // to get VOF DoFHandler and solution

  // Parameters
  const SimulationParameters<dim> &simulation_parameters;

  // Core elements
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl> simulation_control;
  DoFHandler<dim>                    dof_handler;
  std::shared_ptr<FESystem<dim>>     fe;

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

  // Assembler for the matrix and rhs
  std::shared_ptr<PhysicsSubequationsAssemblerBase<ScratchDataType>> assembler;

  // Verbosity
  const Parameters::Verbosity linear_solver_verbosity;
  const Parameters::Verbosity subequation_verbosity;
};

#endif
