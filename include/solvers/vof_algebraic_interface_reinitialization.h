// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_algebraic_interface_reinitialization_h
#define lethe_vof_algebraic_interface_reinitialization_h

#include <solvers/vof_linear_subequations_solver.h>

/**
 * @brief VOF algebraic reinitialization solver.
 *
 * @tparam dim Number of dimensions of the problem.
 */
template <int dim>
class VOFAlgebraicInterfaceReinitialization : public VOFLinearSubequationsSolver<dim>
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
   * @param[in] p_multiphysics_interface Multiphysics interface object used to
   * get information from physics.
   *
   * @param[in,out] p_subequations_interface Subequations interface object used
   * to get information from other subequations and store information from the
   * current one.
   */
  VOFAlgebraicInterfaceReinitialization(
    const SimulationParameters<dim> &p_simulation_parameters,
    const ConditionalOStream        &p_pcout,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                  &p_triangulation,
    MultiphysicsInterface<dim>    *p_multiphysics_interface,
    VOFSubequationsInterface<dim> *p_subequations_interface)
    : VOFLinearSubequationsSolver<dim>(
        VOFSubequationsID::curvature_projection,
        p_simulation_parameters,
        p_simulation_parameters.multiphysics.vof_parameters
          .surface_tension_force.verbosity,
        p_pcout,
        p_triangulation,
        p_multiphysics_interface,
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
  ~VOFAlgebraicInterfaceReinitialization() = default;

private:
  /**
   * @brief Assemble system matrix and right-hand side (rhs).
   */
  void
  assemble_system_matrix_and_rhs() override;
};

#endif
