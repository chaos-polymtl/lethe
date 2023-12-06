/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 *
 * Scratch data for the VOF auxiliary physics
 */

#include <core/multiphysics.h>

#include <solvers/multiphysics_interface.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>

#ifndef lethe_VOF_scratch_data_h
#  define lethe_VOF_scratch_data_h

using namespace dealii;


/**
 * @brief Class that stores the information required by the assembly procedure
 * for a VOF free surface equation. Consequently, this class
 * calculates the phase values (values, gradients, laplacians) and the shape
 * method (values, gradients, laplacians) at all the gauss points for all
 * degrees of freedom and stores it into arrays.
 * This class serves as a separation between the evaluation at the gauss point
 * of the variables of interest and their use in the assembly, which is carried
 * out by the assembler methods.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 **/

template <int dim>
class VOFScratchData
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evaluation, since this needs to be done at the cell level.
   *
   * @param properties_manager The physical properties Manager (see physical_properties_manager.h)
   *
   * @param fe_vof The FESystem used to solve the VOF equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  VOFScratchData(const PhysicalPropertiesManager properties_manager,
                 const FiniteElement<dim>       &fe_vof,
                 const Quadrature<dim>          &quadrature,
                 const Mapping<dim>             &mapping,
                 const FiniteElement<dim>       &fe_fd)
    : properties_manager(properties_manager)
    , fe_values_vof(mapping,
                    fe_vof,
                    quadrature,
                    update_values | update_gradients |
                      update_quadrature_points | update_hessians |
                      update_JxW_values)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
  {
    allocate();
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd The scratch data
   */
  VOFScratchData(const VOFScratchData<dim> &sd)
    : fe_values_vof(sd.fe_values_vof.get_mapping(),
                    sd.fe_values_vof.get_fe(),
                    sd.fe_values_vof.get_quadrature(),
                    update_values | update_gradients |
                      update_quadrature_points | update_hessians |
                      update_JxW_values)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values | update_gradients)
  {
    allocate();
  }


  /** @brief Allocates the memory for the scratch
   *
   * This method allocates the necessary memory for all members of the scratch
   *
   */
  void
  allocate();

  /** @brief Reinitialize the content of the scratch
   *
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for the VOF
   *
   * @param previous_solutions The solutions at the previous time steps
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions)
  {
    fe_values_vof.reinit(cell);
    this->quadrature_points = fe_values_vof.get_quadrature_points();
    auto &fe_vof            = fe_values_vof.get_fe();

    if (dim == 2)
      this->cell_size = std::sqrt(4. * cell->measure() / M_PI) / fe_vof.degree;
    else if (dim == 3)
      this->cell_size =
        pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_vof.degree;

    fe_values_vof.get_function_values(current_solution,
                                      this->present_phase_values);
    fe_values_vof.get_function_gradients(current_solution,
                                         this->phase_gradients);
    fe_values_vof.get_function_laplacians(current_solution,
                                          this->phase_laplacians);

    fe_values_vof.get_function_gradients(previous_solutions[0],
                                         this->previous_phase_gradients);


    // Gather previous vof values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        fe_values_vof.get_function_values(previous_solutions[p],
                                          this->previous_phase_values[p]);
      }

    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->JxW[q] = fe_values_vof.JxW(q);

        for (unsigned int k = 0; k < this->n_dofs; ++k)
          {
            // Shape function
            this->phi[q][k]           = fe_values_vof.shape_value(k, q);
            this->grad_phi[q][k]      = fe_values_vof.shape_grad(k, q);
            this->hess_phi[q][k]      = fe_values_vof.shape_hessian(k, q);
            this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
          }
      }
  }

  /** @brief Reinitialize the velocity, calculated by the Fluid Dynamics
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell for which the velocity is reinitialized
   * This cell must be compatible with the Fluid Dynamics FE
   *
   * @param current_solution The present value of the solution for [u,p]
   *
   */

  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType              &current_solution,
                  const std::vector<VectorType> &previous_solutions,
                  const Parameters::ALE<dim>    &ale)
  {
    fe_values_fd.reinit(cell);

    fe_values_fd[velocities_fd].get_function_values(current_solution,
                                                    velocity_values);
    fe_values_fd[velocities_fd].get_function_gradients(
      current_solution, velocity_gradient_values);

    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->velocity_divergences[q] =
          trace(this->velocity_gradient_values[q]);
      }

    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        fe_values_fd[velocities_fd].get_function_values(
          previous_solutions[p], this->previous_velocity_values[p]);
      }

    if (!ale.enabled())
      return;

    // ALE enabled, so extract the ALE velocity and subtract it from the
    // velocity obtained from the fluid dynamics
    Tensor<1, dim>                                  velocity_ale;
    std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
      ale.velocity;
    Vector<double> velocity_ale_vector(dim);

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        velocity_ale_function->vector_value(quadrature_points[q],
                                            velocity_ale_vector);
        for (unsigned int d = 0; d < dim; ++d)
          velocity_ale[d] = velocity_ale_vector[d];

        velocity_values[q] -= velocity_ale;

        for (unsigned int p = 0; p < previous_solutions.size(); ++p)
          {
            this->previous_velocity_values[p][q] -= velocity_ale;
          }
      }
  }

  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;

  // FEValues for the VOF problem
  FEValues<dim> fe_values_vof;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // VOF values
  std::vector<double>         present_phase_values;
  std::vector<Tensor<1, dim>> phase_gradients;
  std::vector<Tensor<1, dim>> previous_phase_gradients;

  std::vector<double>              phase_laplacians;
  std::vector<std::vector<double>> previous_phase_values;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi;
  std::vector<std::vector<double>>         laplacian_phi;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValues<dim> fe_values_fd;

  FEValuesExtractors::Vector velocities_fd;
  // This FEValues must be instantiated for the velocity
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<Tensor<2, dim>>              velocity_gradient_values;
  std::vector<double>                      velocity_divergences;
};

#endif
