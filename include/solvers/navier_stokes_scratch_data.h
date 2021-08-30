/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <core/bdf.h>
#include <core/parameters.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>


#ifndef lethe_navier_stokes_scratch_data_h
#  define lethe_navier_stokes_scratch_data_h

using namespace dealii;


/**
 * @brief NavierStokesScratchData class
 * stores the information required by the assembly procedure
 * for a Navier-Stokes equation. Consequently, this class calculates
 * the velocity (values, gradients, laplacians) and the shape function
 * (values, gradients, laplacians) at all the gauss points for all degrees
 * of freedom and stores it into arrays. Additionnaly, the use can request
 * that this class gathers additional fields for physics which are coupled
 * to the Navier-Stokes equation, such as the free surface. This class
 * serves as a seperation between the evaluation at the gauss point of the
 * variables of interest and their use in the assembly, which is carried out
 * by the assembler functions. For more information on this design, the reader
 * can consult deal.II step-9
 * "https://www.dealii.org/current/doxygen/deal.II/step_9.html". In this latter
 * example, the scratch is a struct instead of a templated class because of the
 * simplicity of step-9.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *  @ingroup solvers
 **/

template <int dim>
class NavierStokesScratchData
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evalution, since this needs to be done at the cell level.
   *
   * @param fe The FESystem used to solve the Navier-Stokes equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes equations are solved
   *
   */
  NavierStokesScratchData(const FESystem<dim> &  fe,
                          const Quadrature<dim> &quadrature,
                          const Mapping<dim> &   mapping)
    : fe_values(mapping,
                fe,
                quadrature,
                update_values | update_quadrature_points | update_JxW_values |
                  update_gradients | update_hessians)
  {
    allocate();

    // By default, the assembly of variables belonging to auxiliary physics is
    // disabled.
    gather_free_surface = false;
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   *  This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param fe The FESystem used to solve the Navier-Stokes equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes equations are solved
   */
  NavierStokesScratchData(const NavierStokesScratchData<dim> &sd)
    : fe_values(sd.fe_values.get_mapping(),
                sd.fe_values.get_fe(),
                sd.fe_values.get_quadrature(),
                update_values | update_quadrature_points | update_JxW_values |
                  update_gradients | update_hessians)
  {
    allocate();
    if (sd.gather_free_surface)
      enable_free_surface(sd.fe_values_free_surface->get_fe(),
                          sd.fe_values_free_surface->get_quadrature(),
                          sd.fe_values_free_surface->get_mapping());
  }


  /** @brief Allocates the memory for the scratch
   *
   * This function allocates the necessary memory for all members of the scratch
   *
   */
  void
  allocate();

  /** @brief Reinitialize the content of the scratch
   *
   * Using the FeValues and the content ofthe solutions, previous solutions and
   * solutions stages, fills all of the class member of the scratch
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for [u,p]
   *
   * @param previous_solutions The solutions at the previous time steps
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods)
   *
   * @param forcing_function The function describing the momentum/mass source term
   *
   * @param beta_force The additional force for flow control. TODO : Deprecate this argument and pass it
   * to the constructor of the assembler
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType &                                    current_solution,
         const std::vector<VectorType> &previous_solutions,
         const std::vector<VectorType> &solution_stages,
         Function<dim> *                forcing_function,
         Tensor<1, dim>                 beta_force)
  {
    this->fe_values.reinit(cell);

    quadrature_points = this->fe_values.get_quadrature_points();
    auto &fe          = this->fe_values.get_fe();

    forcing_function->vector_value_list(quadrature_points, this->rhs_force);

    // Establish the force vector
    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (int d = 0; d < dim; ++d)
          {
            const unsigned int component_i =
              fe.system_to_component_index(d).first;
            this->force[q][d] = this->rhs_force[q](component_i);
          }
        // Correct force to include the dynamic forcing term for flow
        // control
        force[q] = force[q] + beta_force;
      }

    if (dim == 2)
      this->cell_size = std::sqrt(4. * cell->measure() / M_PI) / fe.degree;
    else if (dim == 3)
      this->cell_size = pow(6 * cell->measure() / M_PI, 1. / 3.) / fe.degree;

    // Gather velocity (values, gradient and laplacian)
    this->fe_values[velocities].get_function_values(current_solution,
                                                    this->velocity_values);
    this->fe_values[velocities].get_function_gradients(
      current_solution, this->velocity_gradients);
    this->fe_values[velocities].get_function_laplacians(
      current_solution, this->velocity_laplacians);
    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->velocity_divergences[q] = trace(this->velocity_gradients[q]);
      }

    // Gather previous velocities
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values[velocities].get_function_values(
          previous_solutions[p], previous_velocity_values[p]);
      }

    // Gather velocity stages
    for (unsigned int s = 0; s < solution_stages.size(); ++s)
      {
        this->fe_values[velocities].get_function_values(
          solution_stages[s], stages_velocity_values[s]);
      }


    // Gather pressure (values, gradient)
    fe_values[pressure].get_function_values(current_solution,
                                            this->pressure_values);
    fe_values[pressure].get_function_gradients(current_solution,
                                               this->pressure_gradients);


    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values.JxW(q);
        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Velocity
            this->phi_u[q][k] = this->fe_values[velocities].value(k, q);
            this->div_phi_u[q][k] =
              this->fe_values[velocities].divergence(k, q);
            this->grad_phi_u[q][k] = this->fe_values[velocities].gradient(k, q);
            this->hess_phi_u[q][k] = this->fe_values[velocities].hessian(k, q);
            for (int d = 0; d < dim; ++d)
              this->laplacian_phi_u[q][k][d] = trace(this->hess_phi_u[q][k][d]);
            // Pressure
            this->phi_p[q][k]      = this->fe_values[pressure].value(k, q);
            this->grad_phi_p[q][k] = this->fe_values[pressure].gradient(k, q);
          }
      }
  }


  /**
   * @brief enable_free_surface Enables the collection of the free surface data by the scratch
   *
   * @param fe FiniteElement associated with the free surface.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */

  void
  enable_free_surface(const FiniteElement<dim> &fe,
                      const Quadrature<dim> &   quadrature,
                      const Mapping<dim> &      mapping)
  {
    gather_free_surface    = true;
    fe_values_free_surface = std::make_shared<FEValues<dim>>(
      mapping, fe, quadrature, update_values | update_gradients);

    // Free surface
    phase_values = std::vector<double>(this->n_q_points);
    previous_phase_values =
      std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                       std::vector<double>(this->n_q_points));
    phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  }

  /** @brief Reinitialize the content of the scratch for the free surface
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the free surface FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [alpha]
   *
   * @param previous_solutions The solutions at the previous time steps for [alpha]
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods) for [alpha]
   *
   */

  template <typename VectorType>
  void
  reinit_free_surface(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType &                                    current_solution,
    const std::vector<VectorType> &                       previous_solutions,
    const std::vector<VectorType> & /*solution_stages*/)
  {
    this->fe_values_free_surface->reinit(cell);
    // Gather phase fraction (values, gradient)
    this->fe_values_free_surface->get_function_values(current_solution,
                                                      this->phase_values);
    this->fe_values_free_surface->get_function_gradients(
      current_solution, this->phase_gradient_values);

    // Gather previous phase fraction values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_free_surface->get_function_values(
          previous_solutions[p], previous_phase_values[p]);
      }
  }

  // FEValues for the Navier-Stokes problem
  FEValues<dim>              fe_values;
  unsigned int               n_dofs;
  unsigned int               n_q_points;
  double                     cell_size;
  FEValuesExtractors::Vector velocities;
  FEValuesExtractors::Scalar pressure;

  std::vector<Vector<double>> rhs_force;
  Tensor<1, dim>              beta_force;
  std::vector<Tensor<1, dim>> force;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Velocity and pressure values
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<double>                      velocity_divergences;
  std::vector<Tensor<2, dim>>              velocity_gradients;
  std::vector<Tensor<1, dim>>              velocity_laplacians;
  std::vector<double>                      pressure_values;
  std::vector<Tensor<1, dim>>              pressure_gradients;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> stages_velocity_values;


  // Shape functions
  std::vector<std::vector<double>>         div_phi_u;
  std::vector<std::vector<Tensor<1, dim>>> phi_u;
  std::vector<std::vector<Tensor<3, dim>>> hess_phi_u;
  std::vector<std::vector<Tensor<1, dim>>> laplacian_phi_u;
  std::vector<std::vector<Tensor<2, dim>>> grad_phi_u;
  std::vector<std::vector<double>>         phi_p;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_p;


  /**
   * Scratch component for the free surface auxiliary physics
   */
  bool                             gather_free_surface;
  unsigned int                     n_dofs_free_surface;
  std::vector<double>              phase_values;
  std::vector<std::vector<double>> previous_phase_values;
  std::vector<Tensor<1, dim>>      phase_gradient_values;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_free_surface;
};

#endif
