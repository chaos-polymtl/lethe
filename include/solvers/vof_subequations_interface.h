// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_subequations_interface_h
#define lethe_subequations_interface_h

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_subequations_solver.h>
#include <solvers/vof_assemblers.h>
#include <solvers/vof_scratch_data.h>
#include <solvers/vof_subequations.h>

using namespace dealii;

DeclExceptionMsg(
  ValidityMapNotReset,
  "The subequation validity map is not rest. Call "
  "'set_vof_dof_handler_and_filtered_solution' before 'set_vof_solution'.");

DeclException1(
  InvalidSubequationSolution,
  std::string,
  "The current solution of the "
    << arg1
    << "is invalid. \n"
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
                                             &p_triangulation,
    const std::shared_ptr<SimulationControl> &p_simulation_control);

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
   * @brief Get a pointer to the DoFHandler of the VOF auxiliary physics.
   *
   * @return Pointer to the DoFHandler of the VOF auxiliary physics associated
   * to the set Filtered VOF solution field.
   */
  DoFHandler<dim> *
  get_vof_dof_handler()
  {
    return this->dof_handler_vof;
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
    AssertThrow(this->get_solution_validity(subequation_id),
                InvalidSubequationSolution(
                  this->get_subequation_string(subequation_id)));

    return this->subequations_solutions[subequation_id];
  }

  /**
   * @brief Get a pointer to the set filtered phase fraction solution vector.
   *
   * @return Pointer to the set filtered phase fraction solution vector.
   */
  GlobalVectorType *
  get_vof_filtered_solution()
  {
    return &this->vof_filtered_solution_vector;
  }

  /**
   * @brief Get a pointer to the set phase fraction solution vector.
   *
   * @return Pointer to the set phase fraction solution vector.
   */
  GlobalVectorType *
  get_vof_solution()
  {
    return &this->vof_solution_vector;
  }

  /**
   * @brief Check and return a boolean indicating if the solution of the
   * subequation corresponds to the set @p vof_filtered_solution_vector.
   *
   * @param[in] subequation_id Identifier associated with a specific
   * subequation.
   *
   * @return @p true if the solution of the subequation corresponds to the set
   * @p vof_filtered_solution_vector. Otherwise @p false.
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

  /**
   * @brief Set a new VOF filtered solution vector and DoFHandler for the
   * subequations to be solved and reset their validity.
   *
   * @param[in] vof_filtered_solution_vector Filtered VOF solution vector
   * associated with the subequations to be solved.
   *
   * @paramp[in] dof_handler_vof Pointer to the DoFHandler associated with the
   * specified VOF solution.
   */
  void
  set_vof_filtered_solution_and_dof_handler(
    const GlobalVectorType &vof_filtered_solution_vector,
    DoFHandler<dim>        *dof_handler_vof)
  {
    // Clear all and reset validity map
    this->vof_filtered_solution_vector.clear();
    this->vof_solution_vector.clear();
    reset_subequations_solutions_validity();

    // Set solution values and DoFHandler
    this->vof_filtered_solution_vector = vof_filtered_solution_vector;
    this->dof_handler_vof              = dof_handler_vof;
  }

  /**
   * @brief Set a new VOF solution vector and DoFHandler for the subequations to
   * be solved.
   *
   * @param[in] vof_solution_vector VOF solution vector associated with the
   * subequations to be solved.
   *
   * @note At the moment, this solution is only used in
   * VOFAlgebraicInterfaceReinitialization.
   */
  void
  set_vof_solution(const GlobalVectorType &vof_solution_vector)
  {
    // By checking if the validity map has been reset, we make sure that the VOF
    // solution is consistent with the set filtered solution and the DoFHandler.
    AssertThrow(validity_map_has_been_reset(), ValidityMapNotReset());
    this->vof_solution_vector = vof_solution_vector;
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
   * @brief Check if subequations_solutions_validity has been reset.
   *
   * @return @p true if all booleans in the map are set to @p false and @p false
   * if at least one of them is set to @p true
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
   * subequations to @p false (default) or @p true if specified.
   *
   * @param[in] reset_value Either @p true or @false (default).
   *
   * @note At the moment, reset_value is only set to @p true for the purpose of
   * outputting the initial condition when user requests it.
   */
  void
  reset_subequations_solutions_validity(const bool reset_value = false)
  {
    for (const VOFSubequationsID &subequation_id : this->active_subequations)
      {
        subequations_solutions_validity[subequation_id] = reset_value;
      }
  }

  // Parallel consol output object
  const ConditionalOStream pcout;

  // VOF DoFHandler associated with solved equations
  DoFHandler<dim> *dof_handler_vof;

  // Filtered VOF solution field associated with solved equations
  GlobalVectorType vof_filtered_solution_vector;

  /** VOF solution field associated with solved equations
   * @note Used in algebraic reinitialization for debuging puposes */
  GlobalVectorType vof_solution_vector;

  // Data structure that stores all enabled subequations
  std::vector<VOFSubequationsID> active_subequations;

  /** Subequations are stored within a map of shared pointer to ensure proper
   * memory management */
  std::map<VOFSubequationsID, std::shared_ptr<PhysicsSubequationsSolverBase>>
    subequations;

  /// Map of DoFHandler of subequations
  std::map<VOFSubequationsID, DoFHandler<dim> *> subequations_dof_handler;

  /// Present solutions map
  std::map<VOFSubequationsID, GlobalVectorType *> subequations_solutions;

  /** Booleans indicating if the subequation has been solved with the set VOF
   * solution vector */
  std::map<VOFSubequationsID, bool> subequations_solutions_validity;
};

#endif
