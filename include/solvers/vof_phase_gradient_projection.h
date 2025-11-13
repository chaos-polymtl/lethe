// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_phase_gradient_projection_h
#define lethe_vof_phase_gradient_projection_h

#include <solvers/vof_linear_subequations_solver.h>

#include <deal.II/fe/fe_simplex_p.h>

DeclException1(
  SameFilteredVOFSolution,
  std::string,
  "A new VOF filtered phase fraction solution has not been set. There is no need "
  "to solve once more the following equation: "
    << arg1
    << ". If you wish to solve the subequation for a new filtered phase "
       "fraction field, please set a new VOF filtered solution field with "
       "VOFSubequationsInterface<dim>::set_vof_filtered_solution_and_dof_handler() "
       "before solving the subequation.");

DeclException1(
  NoFilteredVOFSolution,
  std::string,
  "No VOF filtered phase fraction solution has been set. A valid VOF filtered "
  "phase fraction solution is required to solve the "
    << arg1
    << ". Please set a new VOF filtered phase fraction solution field with "
       "VOFSubequationsInterface<dim>::set_vof_filtered_solution_and_dof_handler() "
       "before solving the subequation.");

/**
 * @brief VOF phase fraction gradient L2 projection solver.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class VOFPhaseGradientProjection : public VOFLinearSubequationsSolver<dim>
{
public:
  /**
   * @brief Constructor for the L2 projection of the VOF phase fraction gradient
   * (pfg).
   *
   * @param[in] p_simulation_parameters Simulation parameters.
   *
   * @param[in] p_pcout Parallel cout used to print the information.
   *
   * @param[in] p_triangulation Distributed mesh information.
   *
   * @param[in,out] p_subequations_interface Subequations interface object used
   * to get information from other subequations and store information from the
   * current one.
   */
  VOFPhaseGradientProjection(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                   p_triangulation,
    VOFSubequationsInterface<dim> &p_subequations_interface)
    : VOFLinearSubequationsSolver<dim>(
        VOFSubequationsID::phase_gradient_projection,
        p_simulation_parameters,
        ((p_simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.verbosity != Parameters::Verbosity::quiet) ||
         ((p_simulation_parameters.multiphysics.vof_parameters
             .regularization_method.algebraic_interface_reinitialization
             .enable) &&
          (p_simulation_parameters.multiphysics.vof_parameters
             .regularization_method.verbosity !=
           Parameters::Verbosity::quiet))) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet, // Set to verbose if surface tension
                                        // verbosity is enabled or if algebraic
                                        // interface reinitialization is enabled
                                        // and set to verbose
        p_pcout,
        p_triangulation,
        p_subequations_interface)
  {
    if (this->simulation_parameters.mesh.simplex)
      {
        // For simplex meshes
        const FE_SimplexP<dim> subequation_fe(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->fe      = std::make_shared<FESystem<dim>>(subequation_fe, dim);
        this->mapping = std::make_shared<MappingFE<dim>>(*this->fe);

        this->cell_quadrature =
          std::make_shared<QGaussSimplex<dim>>(this->fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        const FE_Q<dim> subequation_fe(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->fe      = std::make_shared<FESystem<dim>>(subequation_fe, dim);
        this->mapping = std::make_shared<MappingQ<dim>>(this->fe->degree);

        this->cell_quadrature =
          std::make_shared<QGauss<dim>>(this->fe->degree + 1);
      }
  }

  /**
   * @brief Default destructor.
   */
  ~VOFPhaseGradientProjection() = default;

private:
  /**
   * @brief Assemble system matrix and right-hand side (rhs).
   */
  void
  assemble_system_matrix_and_rhs() override;

  /**
   * @brief Check if a new VOF filtered solution has been set.
   */
  void
  check_dependencies_validity() override;

  /**
   * @brief For vector values, compute the normalized vector field.
   */
  void
  compute_normalized_vector_solution() override;
};

#endif
