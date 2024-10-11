// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_phase_gradient_projection_h
#define lethe_vof_phase_gradient_projection_h

#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_subequations_solver.h>
#include <solvers/subequations_interface.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_filter.h>
#include <solvers/vof_scratch_data.h>

#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include "subequations_interface.h"

template <int dim>
class VOFPhaseGradientProjection
  : public PhysicsLinearSubequationsSolver<dim, GlobalVectorType>
{
public:
  /**
   * @brief Default constructor for the VOFPhaseGradientProjection.
   * TODO AMISHGA
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_triangulation
   *
   * @param[in] p_simulation_control
   */
  VOFPhaseGradientProjection(
    const ConditionalOStream        &p_pcout,
    SubequationsInterface<dim>      *p_subequations,
    MultiphysicsInterface<dim>      *p_multiphysics,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       &p_triangulation,
    std::shared_ptr<SimulationControl> &p_simulation_control)
    : PhysicsLinearSubequationsSolver<dim, GlobalVectorType>(p_pcout)
    //    PhysicsLinearSubequationsSolver<dim, GlobalVectorType>(false,
    //        p_simulation_parameters.non_linear_solver.at(PhysicsID::VOF))
    , computing_timer(p_triangulation->get_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , subequations(p_subequations)
    , multiphysics(p_multiphysics)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
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

    // Change the behavior of the timer for situations when you don't want
    // outputs
    if (this->simulation_parameters.timer.type == Parameters::Timer::Type::none)
      this->computing_timer.disable_output();
  }

  /**
   * @brief Default destructor.
   */
  ~VOFPhaseGradientProjection() = default;

  /**
   * @brief Set up the DofHandler and the degree of freedom associated with
   * the physics.
   */
  void
  setup_dofs() override;

  /**
   * @brief Solve linear system of equation using a strategy appropriate
   * for the partial differential equation.
   */
  void
  solve_linear_system_and_update_solution() override;

  /**
   * @brief Assemble and solve linear system when the equation to solve is
   * linear without using the non-linear solver interface.
   */
  void
  solve() override;

private:
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
    VOFPhaseGradientProjectionScratchData<dim>           &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);


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
    VOFPhaseGradientProjectionScratchData<dim>           &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Copy local cell matrix information to global matrix
   *
   * @param[in] copy_data Stores the results of the assembly over a cell.
   */
  virtual void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsCopyData &copy_data);

  /**
   * @brief Copy local right-hand side (rhs) information to global rhs.
   *
   * @param[in] copy_data Stores the results of the assembly over a cell.
   */
  virtual void
  copy_local_rhs_to_global_rhs(const StabilizedMethodsCopyData &copy_data);


  TimerOutput computing_timer;

  SubequationsInterface<dim> *subequations;
  MultiphysicsInterface<dim> *multiphysics; // to get VOF DoFHandler

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
  std::shared_ptr<VOFAssemblerPhaseGradientProjection<dim>> assembler;

  // TODO AMISHGA check if this is necessary
  /*// Solution transfer classes
  std::shared_ptr<
    parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>
    solution_transfer;
  std::vector<parallel::distributed::SolutionTransfer<dim,
  GlobalVectorType>>
    previous_solutions_transfer;*/
};



#endif
