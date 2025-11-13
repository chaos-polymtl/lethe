// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_curvature_projection_h
#define lethe_vof_curvature_projection_h

#include <solvers/vof_linear_subequations_solver.h>

#include <deal.II/fe/fe_simplex_p.h>

/**
 * @brief VOF curvature L2 projection solver.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class VOFCurvatureProjection : public VOFLinearSubequationsSolver<dim>
{
public:
  /**
   * @brief Constructor for the L2 projection of the VOF interface curvature.
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
  VOFCurvatureProjection(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                   p_triangulation,
    VOFSubequationsInterface<dim> &p_subequations_interface)
    : VOFLinearSubequationsSolver<dim>(
        VOFSubequationsID::curvature_projection,
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
        this->fe = std::make_shared<FE_SimplexP<dim>>(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->mapping = std::make_shared<MappingFE<dim>>(*this->fe);
        this->cell_quadrature =
          std::make_shared<QGaussSimplex<dim>>(this->fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        this->fe = std::make_shared<FE_Q<dim>>(
          this->simulation_parameters.fem_parameters.VOF_order);
        this->mapping = std::make_shared<MappingQ<dim>>(this->fe->degree);
        this->cell_quadrature =
          std::make_shared<QGauss<dim>>(this->fe->degree + 1);
      }
  }

  /**
   * @brief Default destructor.
   */
  ~VOFCurvatureProjection() = default;

private:
  /**
   * @brief Assemble system matrix and right-hand side (rhs).
   */
  void
  assemble_system_matrix_and_rhs() override;

  /**
   * @brief Check if the phase gradient L2 projection has been solved.
   */
  void
  check_dependencies_validity() override;

  /**
   * @brief For vector values, compute the normalized vector field.
   * @warning The method is only implemented for the projected phase gradient.
   */
  void
  compute_normalized_vector_solution() override
  {
    return; // Do nothing
  }
};

#endif
