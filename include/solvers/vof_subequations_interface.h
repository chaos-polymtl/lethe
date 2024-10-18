// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_subequations_interface_h
#define lethe_subequations_interface_h

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>

using namespace dealii;

/**
 * @brief IDs associated to the different subequations solved in Lethe.
 */
enum VOFSubequationsID : unsigned int
{
  /// VOF phase fraction gradient L2 projection
  phase_gradient_projection = 0
};

/**
 * @brief Interface for secondary equations (subequations) solved within
 * the VOF auxiliary physics.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class VOFSubequationsInterface
{
public:
  /**
   * @brief Constructor of an interface for subequations that require the use of
   * a solver within the span of the VOF auxiliary physics' resolution.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_multiphysics Multiphysics interface object used to get
   * information from physics.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_simulation_control Object responsible for the control of
   * steady-state and transient simulations. Contains all the information
   * related to time stepping and the stopping criteria.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   */
  VOFSubequationsInterface(
    const SimulationParameters<dim> &p_simulation_parameters,
    MultiphysicsInterface<dim>      *p_multiphysics,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       &p_triangulation,
    std::shared_ptr<SimulationControl> &p_simulation_control,
    ConditionalOStream                 &p_pcout);

  /**
   * @brief Default destructor.
   */
  ~VOFSubequationsInterface() = default;

  /**
   * @brief Setup the DofHandler and the degree of freedom associated with the
   * subequations.
   */
  void
  setup_dofs()
  {
    for (const auto &subequation : this->subequations)
      {
        subequation.second->setup_dofs();
      }
  };

  /**
   * @brief Setup the DofHandler and the degree of freedom associated with the
   * physics.
   *
   * @param[in] subequation_id Identifier associated with the subequation wished
   * to solve.
   */
  void
  setup_specific_subequation_dofs(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    this->subequations.find(subequation_id)->second->setup_dofs();
  };

  /**
   * @brief Cast the appropriate scratch data object.
   *
   * @param[in] subequation_id Identifier associated with the subequation wished
   * to solve.
   *
   * @param[in] fe_subequation FiniteElement object used for solving the wished
   * subequation.
   *
   * @param[in] quadrature Quadrature rule used for the assembly of the matrix
   * and the right-hand side.
   *
   * @param[in] mapping Mapping of the domain used when solving the VOF
   * equation.
   *
   * @param[in] fe_input FiniteElement object from the VOF auxiliary physics or
   * another subequation that is used as an input for the resolution of the
   * current subequation.
   *
   * @return Shared pointer to the scratch data object of the appropriate
   * subequation.
   */
  std::shared_ptr<PhysicsScratchDataBase>
  scratch_data_cast(const VOFSubequationsID  &subequation_id,
                    const FiniteElement<dim> &fe_subequation,
                    const Quadrature<dim>    &quadrature,
                    const Mapping<dim>       &mapping,
                    const FiniteElement<dim> &fe_input);

  /**
   * @brief Cast the appropriate assembler object.
   *
   * @param[in] subequation_id Identifier associated with the subequation wished
   * to solve.
   *
   * @param[in] vof_parameters VOF auxiliary physics parameter set.
   *
   * @return Shared pointer to the assembler object of the appropriate
   * subequation.
   */
  template <typename ScratchDataType>
  std::shared_ptr<VOFSubequationAssemblerBase<ScratchDataType>>
  assembler_cast(const VOFSubequationsID &subequation_id,
                 const Parameters::VOF   &vof_parameters)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    if (subequation_id == VOFSubequationsID::phase_gradient_projection)
      return std::make_shared<
        VOFAssemblerPhaseGradientProjection<dim, ScratchDataType>>(
        vof_parameters);
    else // At the moment, only one option is possible. This will change with
         // the addition of other subequations to the interface.
      return std::make_shared<
        VOFAssemblerPhaseGradientProjection<dim, ScratchDataType>>(
        vof_parameters);
  }

  /**
   * @brief Call solving method of active subequations.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  void
  solve(const bool &is_post_mesh_adaptation = false)
  {
    for (const auto &subequation : this->subequations)
      {
        subequation.second->solve(is_post_mesh_adaptation);
      }
  }

  /**
   * @brief Call solving method for a specific subequation.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   *
   * @param[in] subequation_id Identifier associated with the subequation wished
   * to solve.
   */
  void
  solve_specific_subequation(const VOFSubequationsID &subequation_id,
                             const bool &is_post_mesh_adaptation = false)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    this->subequations.find(subequation_id)
      ->second->solve(is_post_mesh_adaptation);
  }

  /**
   * @brief Get vector of active subequations.
   *
   * @return Vector of active subequations identifiers.
   */
  std::vector<VOFSubequationsID>
  get_active_subequations()
  {
    return this->active_subequations;
  }

  /**
   * @brief Get a pointer to the DoFHandler of a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return Pointer to the DoFHandler of the specified subequation.
   */
  DoFHandler<dim> *
  get_dof_handler(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    return this->subequations_dof_handler[subequation_id];
  }

  /**
   * @brief Get a pointer to the solution vector of a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return Pointer to the solution vector of the specified subequation.
   */
  GlobalVectorType *
  get_solution(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    return this->subequations_solutions[subequation_id];
  }

  /**
   * @brief Get the string associated with a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return String associated with the specified subequation.
   *
   * @note This is used for printing purposes.
   */
  std::string
  get_subequation_string(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    std::string subequation_string;

    if (subequation_id == VOFSubequationsID::phase_gradient_projection)
      subequation_string = "VOF phase fraction gradient L2 projection";

    return subequation_string;
  }

  /**
   * @brief Set the DoFHandler associated with a specific subequation in the
   * interface.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @param[in] dof_handler Pointer to the DoFHandler of a specific subequation.
   */
  void
  set_dof_handler(const VOFSubequationsID &subequation_id,
                  DoFHandler<dim>         *dof_handler)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    this->subequations_dof_handler[subequation_id] = dof_handler;
  }

  /**
   * @brief Set the solution associated with a specific subequation in the
   * interface.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @param[in] solution_vector Pointer to the solution vector of a specific
   * subequation.
   */
  void
  set_solution(const VOFSubequationsID &subequation_id,
               GlobalVectorType        *solution_vector)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    this->subequations_solutions[subequation_id] = solution_vector;
  }

private:
  MultiphysicsInterface<dim> *multiphysics;

  ConditionalOStream pcout;

  // Data structure that stores all enabled physics
  std::vector<VOFSubequationsID> active_subequations;

  // Subequations stored within a map of shared pointer to ensure proper
  // memory management.
  std::map<VOFSubequationsID, std::shared_ptr<PhysicsSubequationsSolverBase>>
    subequations;

  std::map<VOFSubequationsID, DoFHandler<dim> *> subequations_dof_handler;

  // Present solutions
  std::map<VOFSubequationsID, GlobalVectorType *> subequations_solutions;
};

#endif
