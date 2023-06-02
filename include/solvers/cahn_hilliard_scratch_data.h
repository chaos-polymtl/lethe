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
 * Scratch data for the cahn hilliard auxiliary physics
 */

#include <core/multiphysics.h>

#include <solvers/physical_properties_manager.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>


#ifndef lethe_cahn_hilliard_scratch_data_h
#  define lethe_cahn_hilliard_scratch_data_h

using namespace dealii;


/**
 * @brief CahnHilliardScratchData class
 * stores the information required by the assembly procedure
 * for the Cahn-Hilliard equations Consequently, this class
 * calculates the phase field parameter Phi (values, gradients, laplacians),
 * the chemical potential eta (values,gradients,laplacians) and the shape function
 * (values, gradients, laplacians) at all the gauss points for all degrees
 * of freedom and stores it into arrays. Additionnaly, the use can request
 * that this class gathers additional fields for physics which are coupled
 * to the Cahn-Hilliard equations, such as the velocity which is required. This class
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
class CahnHilliardScratchData
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evalution, since this needs to be done at the cell level.
   *
   * @param fe_ch The FESystem used to solve the Cahn-Hilliard equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Cahn-Hilliard equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  CahnHilliardScratchData(const PhysicalPropertiesManager &properties_manager,
                    const FESystem<dim> &      fe_ch,
                    const Quadrature<dim> &          quadrature,
                    const Mapping<dim> &             mapping,
                    const FiniteElement<dim> &       fe_fd)
    : properties_manager(properties_manager)
    , fe_values_ch(mapping,
                       fe_ch,
                       quadrature,
                       update_values | update_quadrature_points |
                         update_JxW_values | update_gradients | update_hessians)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values)
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
  CahnHilliardScratchData(const CahnHilliardScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_ch(sd.fe_values_ch.get_mapping(),
                       sd.fe_values_ch.get_fe(),
                       sd.fe_values_ch.get_quadrature(),
                       update_values | update_quadrature_points |
                         update_JxW_values | update_gradients | update_hessians)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values)
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
   * Using the FeValues and the content of the solutions, previous solutions and
   * solutions stages, fills all of the class member of the scratch
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for [Phi,eta]
   *
   * @param previous_solutions The solutions at the previous time steps
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods)
   *
   * @param source_function The function describing the source term in Cahn-Hilliard
   * equations
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType &                                    current_solution,
         const std::vector<VectorType> &previous_solutions,
         const std::vector<VectorType> &solution_stages,
         Function<dim> *                source_function)
  {
    this->fe_values_ch.reinit(cell);

    quadrature_points = this->fe_values_ch.get_quadrature_points();
    auto &fe_ch   = this->fe_values_ch.get_fe();

    source_function->value_list(quadrature_points, source);

    if (dim == 2)
      this->cell_size =
        std::sqrt(4. * cell->measure() / M_PI) / fe_ch.degree;
    else if (dim == 3)
      this->cell_size =
        pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_ch.degree;

    // Gather Phi and eta (values, gradient and laplacian)
    this->fe_values_ch[phase_order].get_function_values(current_solution,
                                               this->phase_order_values);
    this->fe_values_ch[phase_order].get_function_gradients(current_solution,
                                                  this->phase_order_gradients);
    this->fe_values_ch[phase_order].get_function_laplacians(current_solution,
                                                   this->phase_order_laplacians);


    // Gather previous phase order values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_ch[phase_order].get_function_values(previous_solutions[p],
                                                   previous_phase_order_values[p]);
      }

    // Gather phase order stages
    for (unsigned int s = 0; s < solution_stages.size(); ++s)
      {
        this->fe_values_ch[phase_order].get_function_values(solution_stages[s],
                                                   stages_phase_order_values[s]);
      }

    // Gather previous chemical potential values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_ch[chemical_potential].get_function_values(previous_solutions[p],
                                                   previous_chemical_potential_values[p]);
      }

    // Gather chemical potential stages
    for (unsigned int s = 0; s < solution_stages.size(); ++s)
      {
        this->fe_values_ch[chemical_potential].get_function_values(solution_stages[s],
                                                   stages_chemical_potential_values[s]);
      }


    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_ch.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Shape functions for the phase order
            this->phi_phase[q][k]      = this->fe_values_ch[phase_order].value(k, q);
            this->grad_phi_phase[q][k] = this->fe_values_ch[phase_order].gradient(k, q);
            this->hess_phi_phase[q][k] = this->fe_values_ch[phase_order].hessian(k, q);
            this->laplacian_phi_phase[q][k] = trace(this->hess_phi_phase[q][k]);

            // Shape functions for the chemical potential

            this->phi_potential[q][k]      = this->fe_values_ch[chemical_potential].value(k, q);
            this->grad_phi_potential[q][k] = this->fe_values_ch[chemical_potential].gradient(k, q);
            this->hess_phi_potential[q][k] = this->fe_values_ch[chemical_potential].hessian(k, q);
            this->laplacian_phi_potential[q][k] = trace(this->hess_phi_potential[q][k]);

          }
      }
  }

  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType &current_solution)
  {
    this->fe_values_fd.reinit(cell);

    this->fe_values_fd[velocities].get_function_values(current_solution,
                                                       velocity_values);
  }

  /** @brief Calculates the physical properties. This function calculates the physical properties
   * that may be required by the Cahn-Hilliard equations. Namely W (what is it?),
   * the mobility function M, the mobility factor D and the interface thickness
   * epsilon.
   *
   */
//  void
//  calculate_physical_properties();

  // Physical properties //Mettre les propriétés de Cahn-Hilliard (W,M,D,epsilon)
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
//  std::vector<double>                  tracer_diffusivity;
//  std::vector<double>                  tracer_diffusivity_0;
//  std::vector<double>                  tracer_diffusivity_1;







  // FEValues for the Cahn-Hilliard problem
  FEValues<dim> fe_values_ch;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;
  FEValuesExtractors::Scalar phase_order;
  FEValuesExtractors::Scalar chemical_potential;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Phase order and chemical potential values
  std::vector<double>              phase_order_values;
  std::vector<Tensor<1, dim>>      phase_order_gradients;
  std::vector<double>              phase_order_laplacians;
  std::vector<std::vector<double>> previous_phase_order_values;
  std::vector<std::vector<double>> stages_phase_order_values;

  std::vector<double>              chemical_potential_values;
  std::vector<Tensor<1, dim>>      chemical_potential_gradients;
  std::vector<double>              chemical_potential_laplacians;
  std::vector<std::vector<double>> previous_chemical_potential_values;
  std::vector<std::vector<double>> stages_chemical_potential_values;

  // Source term
  std::vector<double> source;

  // Shape functions for the phase order and the chemical potential
  std::vector<std::vector<double>>         phi_phase;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi_phase;
  std::vector<std::vector<double>>         laplacian_phi_phase;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_phase;
  std::vector<std::vector<double>>         phi_potential;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi_potential;
  std::vector<std::vector<double>>         laplacian_phi_potential;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_potential;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
};

#endif
