/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_vof_algebraic_reinitialization_h
#define lethe_vof_algebraic_reinitialization_h

#include <solvers/physics_subequations_solver.h>
#include <solvers/vof.h>
#include <solvers/vof_algebraic_reinitialization_assemblers.h>
#include <solvers/vof_algebraic_reinitialization_scratch_data.h>

template <int dim>
class VOFAlgebraicReinitialization
  : public PhysicsSubequationsSolver<dim, GlobalVectorType>
{
public:
  VOFAlgebraicReinitialization(
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control)
    : PhysicsSubequationsSolver<dim, GlobalVectorType>(
        p_simulation_parameters.non_linear_solver.at(PhysicsID::VOF))
    , computing_timer(p_triangulation->get_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , simulation_control(p_simulation_control)
    , dof_handler(*triangulation)
  {
    // TODO AMISHGA check if we keep VOF order
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        fe = std::make_shared<FE_SimplexP<dim>>(
          simulation_parameters.fem_parameters.VOF_order);
        mapping         = std::make_shared<MappingFE<dim>>(*fe);
        cell_quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
        face_quadrature =
          std::make_shared<QGaussSimplex<dim - 1>>(fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe = std::make_shared<FE_Q<dim>>(
          simulation_parameters.fem_parameters.VOF_order);
        mapping         = std::make_shared<MappingQ<dim>>(fe->degree);
        cell_quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
        face_quadrature = std::make_shared<QGauss<dim - 1>>(fe->degree + 1);
      }

    // TODO AMISHGA check if this is necessary
    //  (with bdf might need to store previous solutions when checkpointing and
    //  all)
    /*// Allocate solution transfer
    solution_transfer = std::make_shared<
      parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>(
      dof_handler);

    // Set size of previous solutions using BDF schemes information
    previous_solutions.resize(maximum_number_of_previous_solutions());

    // Prepare previous solutions transfer
    previous_solutions_transfer.reserve(previous_solutions.size());
    for (unsigned int i = 0; i < previous_solutions.size(); ++i)
      {
        previous_solutions_transfer.emplace_back(
          parallel::distributed::SolutionTransfer<dim, GlobalVectorType>(
            this->dof_handler));
      }*/

    // Change the behavior of the timer for situations when you don't want
    // outputs
    if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
      this->computing_timer.disable_output();
  }

  /**
   * @brief Default destructor
   */
  ~VOFAlgebraicReinitialization() = default;


private:
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
   *
   * @param[in] cell
   *
   * @param[in,out] scratch_data
   *
   * @param[in,out] copy_data
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    VOFAlgebraicReinitializationScratchData<dim>         &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   *
   * @param[in] cell
   *
   * @param[in,out] scratch_data
   *
   * @param[in,out] copy_data
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    VOFAlgebraicReinitializationScratchData<dim>         &scratch_data,
    StabilizedMethodsCopyData                            &copy_data);

  /**
   * @brief Sets-up the vector of assembler functions
   */
  virtual void
  setup_assemblers();



  GlobalVectorType nodal_phase_fraction_owned;

  TimerOutput computing_timer; // TODO AMISHGA needed?

  const SimulationParameters<dim> &simulation_parameters;

  // Core elements
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  std::shared_ptr<SimulationControl>  simulation_control;
  DoFHandler<dim>                     dof_handler;
  std::shared_ptr<FiniteElement<dim>> fe;
  //  ConvergenceTable                    error_table; TODO AMISHGA needed?

  // Mapping and Quadrature
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  // Previous solutions vectors
  std::vector<GlobalVectorType> previous_solutions;

  // Solution storage
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  GlobalVectorType evaluation_point;
  GlobalVectorType local_evaluation_point;
  GlobalVectorType newton_update;
  GlobalVectorType present_solution;
  GlobalVectorType system_rhs;
  //  AffineConstraints<double>                          nonzero_constraints;
  //  AffineConstraints<double>                          bounding_constraints;
  //  AffineConstraints<double>                          zero_constraints;
  TrilinosWrappers::SparseMatrix                     system_matrix;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  // TODO AMISHGA check if this is necessary
  /*// Solution transfer classes
  std::shared_ptr<
    parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>
    solution_transfer;
  std::vector<parallel::distributed::SolutionTransfer<dim,
  GlobalVectorType>>
    previous_solutions_transfer;*/


  // Assemblers for the matrix and rhs
  std::vector<std::shared_ptr<VOFAlgebraicReinitializationAssemblerBase<dim>>>
    assemblers;
};


#endif
