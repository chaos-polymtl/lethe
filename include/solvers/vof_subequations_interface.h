// SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_subequations_interface_h
#define lethe_subequations_interface_h

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_subequations.h>

using namespace dealii;

DeclExceptionMsg(
  ValidityMapNotReset,
  "The subequation validity map is not reset. Call "
  "'set_vof_solution_and_dof_handler' before 'set_vof_solution'.");

DeclException1(
  InvalidSubequationSolution,
  std::string,
  "The current solution of the "
    << arg1
    << " is invalid. \n"
       "A new VOF solution has been set but the subequation was not solved "
       "afterwards.");

DeclException2(PhaseGradientProjectionIsInvalid,
               std::string,
               std::string,
               "The " << arg1 << " has to be solved before " << arg2 << ".");

DeclException2(CurvatureProjectionIsInvalid,
               std::string,
               std::string,
               "The " << arg1 << " has to be solved before " << arg2 << ".");

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
   * @param[in] p_pcout Parallel cout used to print the information.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_simulation_control SimulationControl object.
   */
  VOFSubequationsInterface(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                              p_triangulation,
    const std::shared_ptr<SimulationControl> &p_simulation_control)
    : pcout(p_pcout)
  {
    initialize_subequations(p_simulation_parameters,
                            p_triangulation,
                            p_simulation_control);
  }

  /**
   * @brief Default destructor.
   */
  ~VOFSubequationsInterface() = default;

  /**
   * @brief Initialize subequations necessary for the simulation and reset
   * subeqaution validation map.
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in] p_simulation_control SimulationControl object.
   */
  void
  initialize_subequations(
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                              p_triangulation,
    const std::shared_ptr<SimulationControl> &p_simulation_control);

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
   * @brief Call solving method of active subequations.
   */
  void
  solve()
  {
    for (const auto &subequation : this->subequations)
      {
        subequation.second->solve();
      }
  }

  /**
   * @brief Call solving method for a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with the subequation wished
   * to solve.
   */
  void
  solve_specific_subequation(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    this->subequations.find(subequation_id)->second->solve();
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
   * @brief Get a reference to the DoFHandler of a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return Reference to the DoFHandler of the specified subequation.
   */
  const DoFHandler<dim> &
  get_dof_handler(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());

    return *this->subequations_dof_handler[subequation_id];
  }

  /**
   * @brief Get a reference to the DoFHandler of the VOF auxiliary physics.
   *
   * @return Reference to the DoFHandler of the VOF auxiliary physics associated
   * to the set VOF solution field.
   */
  const DoFHandler<dim> &
  get_vof_dof_handler()
  {
    return this->dof_handler_vof->get();
  }

  /**
   * @brief Get a reference to the solution vector of a specific subequation.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return Reference to the solution vector of the specified subequation.
   */
  const GlobalVectorType &
  get_solution(const VOFSubequationsID &subequation_id)
  {
    AssertThrow(this->get_solution_validity(subequation_id),
                InvalidSubequationSolution(
                  this->get_subequation_string(subequation_id)));

    return *this->subequations_solutions[subequation_id];
  }

  /**
   * @brief Get a reference to the phase fraction solution vector.
   *
   * @return Reference to the phase fraction solution vector.
   */
  const GlobalVectorType &
  get_vof_solution()
  {
    return *this->vof_solution_vector;
  }

  /**
   * @brief Check and return a boolean indicating if the solution of the
   * subequation has been solved for the corresponding @p dof_handler_vof
   * and/or @p vof_solution_vector previously set with
   * VOFSubequationsInterface<dim>::set_vof_solution_and_dof_handler
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return
   *  - @p true if the solution of the subequation corresponds to the set
   * @p dof_handler_vof and/or @p vof_solution_vector.
   *  - Otherwise, @p false.
   */
  bool
  get_solution_validity(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    return this->subequations_solutions_validity[subequation_id];
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
    else if (subequation_id == VOFSubequationsID::curvature_projection)
      subequation_string = "VOF curvature L2 projection";
    else if (subequation_id ==
             VOFSubequationsID::algebraic_interface_reinitialization)
      subequation_string = "VOF algebraic interface reinitialization";
    else
      throw(std::invalid_argument("Invalid VOFSubequationID. Options are: \n"
                                  " <phase_gradient_projection>\n"
                                  " <curvature_projection>\n"
                                  " <algebraic_interface_reinitialization>"));

    return subequation_string;
  }

  /**
   * @brief Set the DoFHandler associated with a specific subequation in the
   * interface.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @param[in] dof_handler Shared pointer to the DoFHandler of a specific
   * subequation.
   */
  void
  set_dof_handler(const VOFSubequationsID         &subequation_id,
                  std::shared_ptr<DoFHandler<dim>> dof_handler)
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
   * @param[in] solution_vector Shared pointer to the solution vector of a
   * specific subequation.
   */
  void
  set_solution(const VOFSubequationsID          &subequation_id,
               std::shared_ptr<GlobalVectorType> solution_vector)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    this->subequations_solutions[subequation_id] = solution_vector;
  }

  /**
   * @brief Set a new VOF solution vector and DoFHandler for the subequations
   * to be solved and reset their validity.
   *
   * @param[in] vof_solution_vector VOF solution vector associated with the
   * subequations to be solved.
   *
   * @param[in] dof_handler_vof Reference to the DoFHandler associated with the
   * specified VOF solution.
   */
  void
  set_vof_solution_and_dof_handler(const GlobalVectorType &vof_solution_vector,
                                   const DoFHandler<dim>  &dof_handler_vof)
  {
    // Reset validity map
    reset_subequations_solutions_validity();

    // Set solution values and DoFHandler
    this->vof_solution_vector = std::cref(vof_solution_vector);
    this->dof_handler_vof     = std::cref(dof_handler_vof);
  }

  /**
   * @brief Set the solution as valid after solving it.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   */
  void
  set_solution_valid(const VOFSubequationsID &subequation_id)
  {
    AssertThrow((std::find(this->active_subequations.begin(),
                           this->active_subequations.end(),
                           subequation_id) != this->active_subequations.end()),
                ExcInternalError());
    this->subequations_solutions_validity[subequation_id] = true;
  }

  /**
   * @brief Check if all booleans in subequations_solutions_validity have been
   * reset to @p false.
   *
   * @return
   * - @p true if all booleans in the map are set to @p false or;
   * - @p false if at least one of the booleans is set to @p true
   */
  bool
  validity_map_has_been_reset()
  {
    for (const auto &pair : subequations_solutions_validity)
      {
        if (pair.second == true)
          return false;
      }
    return true;
  }

private:
  /**
   * @brief Reset all boolean in the validity associated with active
   * subequations to @p false.
   */
  void
  reset_subequations_solutions_validity()
  {
    for (const VOFSubequationsID &subequation_id : this->active_subequations)
      {
        subequations_solutions_validity[subequation_id] = false;
      }
  }

  /// Parallel consol output object
  const ConditionalOStream pcout;

  /// VOF DoFHandler associated with solved equations
  std::optional<std::reference_wrapper<const DoFHandler<dim>>> dof_handler_vof;

  /** VOF solution field associated with solved equations
   * @note Serves as initial condition in algebraic reinitialization */
  std::optional<std::reference_wrapper<const GlobalVectorType>>
    vof_solution_vector;

  /// Data structure that stores all enabled subequations
  std::vector<VOFSubequationsID> active_subequations;

  /** Subequations are stored within a map of shared pointer to ensure proper
   * memory management */
  std::map<VOFSubequationsID, std::shared_ptr<PhysicsSubequationsSolverBase>>
    subequations;

  /// Map of DoFHandler of subequations
  std::map<VOFSubequationsID, std::shared_ptr<DoFHandler<dim>>>
    subequations_dof_handler;

  /// Present solutions map
  std::map<VOFSubequationsID, std::shared_ptr<GlobalVectorType>>
    subequations_solutions;

  /**
   * Booleans indicating if the subequation has been solved with the set VOF
   * solution vector. The VOF solution vector is set with
   * VOFSubequationsInterface<dim>::set_vof_solution_and_dof_handler.
   *
   * During initialization
   * (VOFSubequationsInterface<dim>::initialize_subequations),
   * all boolean are set to @p false.
   *
   * We start by setting a VOF solution and DOFHandler
   * (VOFSubequationsInterface<dim>::set_vof_solution_and_dof_handler),
   * and then start solving equations making sure that dependencies are met
   * (VOFPhaseGradientProjection<dim>::check_dependencies_validity,
   * VOFCurvatureProjection<dim>::check_dependencies_validity, and
   * VOFAlgebraicInterfaceReinitialization<dim>::check_dependencies_validity).
   * For instance, to solve the curvature projection subequation, the phase
   * gradient projection subequation has to be solved before as the projected
   * phase gradient is used in the curvature projection subequation.
   *
   * A subequation in the validity map is to @p true only when it is solved, and
   * it remains valid until a new VOF solution and/or DOFHandler are
   * set through
   * VOFSubequationsInterface<dim>::set_vof_solution_and_dof_handler.
   *
   * This validation mechanism allows, to check if:
   * - A new VOF solution has been set;
   * - All dependencies of the subequation have been solved;
   * - The current solution of the subequation is valid, and avoid unnecessarily
   * solving it once more.
   *
   * @note
   * VOFSubequationsInterface<dim>::set_vof_solution_and_dof_handler
   * also resets the validity map values to @p false.
   */
  std::map<VOFSubequationsID, bool> subequations_solutions_validity;
};

#endif
