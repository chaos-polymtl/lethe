// SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_phase_gradient_projection_h
#define lethe_cls_phase_gradient_projection_h

#include <solvers/cls_linear_subequations_solver.h>

#include <deal.II/fe/fe_simplex_p.h>

DeclException1(
  SameCLSSolution,
  std::string,
  "A new CLS phase indicator solution has not been set. There is no need "
  "to solve once more the following equation: "
    << arg1
    << ". If you wish to solve the subequation for a new phase indicator field, "
       "please set a new CLS solution field with "
       "CLSSubequationsInterface<dim>::set_cls_solution_and_dof_handler() "
       "before solving the subequation.");

DeclException1(
  NoCLSSolution,
  std::string,
  "No CLS phase indicator solution has been set. A valid CLS phase indicator "
  "solution is required to solve the "
    << arg1
    << ". Please set a new CLS phase indicator solution field with "
       "CLSSubequationsInterface<dim>::set_cls_solution_and_dof_handler() "
       "before solving the subequation.");

/**
 * @brief CLS phase indicator gradient L2 projection solver.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class CLSPhaseGradientProjection : public CLSLinearSubequationsSolver<dim>
{
public:
  /**
   * @brief Constructor for the L2 projection of the CLS phase indicator gradient
   * (PIG).
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
  CLSPhaseGradientProjection(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                   p_triangulation,
    CLSSubequationsInterface<dim> &p_subequations_interface)
    : CLSLinearSubequationsSolver<dim>(
        CLSSubequationsID::phase_gradient_projection,
        p_simulation_parameters,
        ((p_simulation_parameters.multiphysics.cls_parameters
            .surface_tension_force.verbosity != Parameters::Verbosity::quiet) ||
         ((p_simulation_parameters.multiphysics.cls_parameters
             .reinitialization_method.pde_based_interface_reinitialization
             .enable) &&
          (p_simulation_parameters.multiphysics.cls_parameters
             .reinitialization_method.verbosity !=
           Parameters::Verbosity::quiet))) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet, // Set to verbose if surface tension
                                        // verbosity is enabled or if PDE-based
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
          this->simulation_parameters.fem_parameters.CLS_order);
        this->fe      = std::make_shared<FESystem<dim>>(subequation_fe, dim);
        this->mapping = std::make_shared<MappingFE<dim>>(*this->fe);

        this->cell_quadrature =
          std::make_shared<QGaussSimplex<dim>>(this->fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        const FE_Q<dim> subequation_fe(
          this->simulation_parameters.fem_parameters.CLS_order);
        this->fe      = std::make_shared<FESystem<dim>>(subequation_fe, dim);
        this->mapping = std::make_shared<MappingQ<dim>>(this->fe->degree);

        this->cell_quadrature =
          std::make_shared<QGauss<dim>>(this->fe->degree + 1);
      }
  }

  /**
   * @brief Default destructor.
   */
  ~CLSPhaseGradientProjection() = default;

private:
  /**
   * @brief Assemble system matrix and right-hand side (rhs).
   */
  void
  assemble_system_matrix_and_rhs() override;

  /**
   * @brief Check if a new CLS solution has been set.
   */
  void
  check_dependencies_validity() override;
};

#endif
