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
 * @brief VOFScratchData class
 * stores the information required by the assembly procedure
 * for a VOF free surface equation. Consequently, this class
 * calculates the phase values (values, gradients, laplacians) and the shape
 * function (values, gradients, laplacians) at all the gauss points for all
 * degrees of freedom and stores it into arrays.
 * This class serves as a seperation between the evaluation at the gauss point
 *of the variables of interest and their use in the assembly, which is carried
 *out by the assembler functions.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *  @ingroup solvers
 **/

template <int dim>
class VOFScratchData
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
  VOFScratchData(const FiniteElement<dim> &fe_fs,
                 const Quadrature<dim> &   quadrature,
                 const Mapping<dim> &      mapping,
                 const FiniteElement<dim> &fe_navier_stokes)
    : fe_values_fs(mapping,
                   fe_fs,
                   quadrature,
                   update_values | update_gradients | update_quadrature_points |
                     update_hessians | update_JxW_values)
    , fe_values_navier_stokes(mapping,
                              fe_navier_stokes,
                              quadrature,
                              update_values | update_gradients)
  {
    allocate();
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
  VOFScratchData(const VOFScratchData<dim> &sd)
    : fe_values_fs(sd.fe_values_fs.get_mapping(),
                   sd.fe_values_fs.get_fe(),
                   sd.fe_values_fs.get_quadrature(),
                   update_values | update_gradients | update_quadrature_points |
                     update_hessians | update_JxW_values)
    , fe_values_navier_stokes(sd.fe_values_navier_stokes.get_mapping(),
                              sd.fe_values_navier_stokes.get_fe(),
                              sd.fe_values_navier_stokes.get_quadrature(),
                              update_values | update_gradients)
  {
    allocate();
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
   * @param source_function The function describing the VOF source term
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType &                                    current_solution,
         const std::vector<VectorType> &previous_solutions,
         const std::vector<VectorType> &solution_stages)
  {
    this->fe_values_fs.reinit(cell);
    quadrature_points = this->fe_values_fs.get_quadrature_points();
    auto &fe_fs       = this->fe_values_fs.get_fe();

    if (dim == 2)
      this->cell_size = std::sqrt(4. * cell->measure() / M_PI) / fe_fs.degree;
    else if (dim == 3)
      this->cell_size = pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_fs.degree;

    this->fe_values_fs.get_function_values(current_solution,
                                           this->present_phase_values);
    this->fe_values_fs.get_function_gradients(current_solution,
                                              this->phase_gradients);
    this->fe_values_fs.get_function_laplacians(current_solution,
                                               this->phase_laplacians);


    // Gather previous fs values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_fs.get_function_values(previous_solutions[p],
                                               previous_phase_values[p]);
      }

    // Gather fs stages
    for (unsigned int s = 0; s < solution_stages.size(); ++s)
      {
        this->fe_values_fs.get_function_values(solution_stages[s],
                                               stages_phase_values[s]);
      }


    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_fs.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Shape function
            this->phi[q][k]           = this->fe_values_fs.shape_value(k, q);
            this->grad_phi[q][k]      = this->fe_values_fs.shape_grad(k, q);
            this->hess_phi[q][k]      = this->fe_values_fs.shape_hessian(k, q);
            this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
          }
      }
  }

  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType &current_solution)
  {
    this->fe_values_navier_stokes.reinit(cell);

    this->fe_values_navier_stokes[velocities].get_function_values(
      current_solution, velocity_values);
    this->fe_values_navier_stokes[velocities].get_function_gradients(
      current_solution, velocity_gradient_values);
  }

  // FEValues for the VOF problem
  FEValues<dim> fe_values_fs;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Free surface values
  std::vector<double>              present_phase_values;
  std::vector<Tensor<1, dim>>      phase_gradients;
  std::vector<double>              phase_laplacians;
  std::vector<std::vector<double>> previous_phase_values;
  std::vector<std::vector<double>> stages_phase_values;

  // double velocity_fem_degree =
  //  simulation_parameters.fem_parameters.velocity_order;

  // Source term
  std::vector<double> source;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi;
  std::vector<std::vector<double>>         laplacian_phi;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValues<dim> fe_values_navier_stokes;

  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<Tensor<2, dim>> velocity_gradient_values;
};

#endif
