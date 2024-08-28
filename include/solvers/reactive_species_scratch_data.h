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
 * Scratch data for the ReactiveSpecies auxiliary physics
 */


#ifndef lethe_reactive_species_scratch_data_h
#define lethe_reactive_species_scratch_data_h

#include <core/multiphysics.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physical_properties_manager.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>


using namespace dealii;


/**
 * @brief Class that stores the information required by the assembly procedure.....
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 **/

template <int dim>
class ReactiveSpeciesScratchData
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evaluation, since this needs to be done at the cell level.
   *
   * @param properties_manager The physical properties manager (see physical_properties_manager.h)
   *
   * @param fe_reactive_species The FESystem used to solve the Reactive species equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Reactive species equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  ReactiveSpeciesScratchData(
    const PhysicalPropertiesManager &properties_manager,
    const FESystem<dim>             &fe_reactive_species,
    const Quadrature<dim>           &quadrature,
    const Mapping<dim>              &mapping,
    const FiniteElement<dim>        &fe_fd)
    : properties_manager(properties_manager)
    , fe_values_extractors(4) // TODO Make flexible with number of species
    , fe_values_reactive_species(mapping,
                                 fe_reactive_species,
                                 quadrature,
                                 update_values | update_quadrature_points |
                                   update_JxW_values | update_gradients |
                                   update_hessians)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values)
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
   * @param fe The FESystem used to solve the Navier-Stokes equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes equations are solved
   */
  ReactiveSpeciesScratchData(const ReactiveSpeciesScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_extractors() // TODO Make flexible with number of species
    , fe_values_reactive_species(sd.fe_values_reactive_species.get_mapping(),
                                 sd.fe_values_reactive_species.get_fe(),
                                 sd.fe_values_reactive_species.get_quadrature(),
                                 update_values | update_quadrature_points |
                                   update_JxW_values | update_gradients |
                                   update_hessians)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values)
  {
    // TODO Make flexible with number of species
    for (unsigned int i = 0; i < 4; i++)
      {
        fe_values_extractors.emplace_back(i);
      }
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
   * Using the FeValues and the content of the solutions and previous solutions
   * , fills all of the class member of the scratch
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for [Phi,eta]
   *
   * @param previous_solutions The solutions at the previous time steps
   *
   * @param source_function The function describing the source term in Reactive species
   * equations
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         Function<dim>                 *source_function)
  {
    // TODO Flexible number of species
    for (unsigned int i = 0; i < 4; i++)
      {
        this->fe_values_extractors[i].component = 0;
      }

    this->fe_values_reactive_species.reinit(cell);

    quadrature_points =
      this->fe_values_reactive_species.get_quadrature_points();
    auto &fe_reactive_species = this->fe_values_reactive_species.get_fe();

    // TODO Flexible number of species
    for (unsigned int i = 0; i < 4; i++)
      {
        source_function->value_list(quadrature_points, source[i], i);
      }

    // TODO Use the new utilities function for that
    if (dim == 2)
      this->cell_size =
        std::sqrt(4. * cell->measure() / M_PI) / fe_reactive_species.degree;
    else if (dim == 3)
      this->cell_size =
        pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_reactive_species.degree;

    // Gather Phi and eta (values, gradient and laplacian)
    for (unsigned int i = 0; i < 4; i++)
      {
        this->fe_values_reactive_species[fe_values_extractors[i]]
          .get_function_values(current_solution,
                               this->reactive_species_values[i]);
        this->fe_values_reactive_species[fe_values_extractors[i]]
          .get_function_gradients(current_solution,
                                  this->reactive_species_gradients[i]);
        this->fe_values_reactive_species[fe_values_extractors[i]]
          .get_function_laplacians(current_solution,
                                   this->reactive_species_laplacians[i]);
      }


    // Gather previous phase order values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        for (unsigned int i = 0; i < 4; i++)
          {
            this->fe_values_reactive_species[fe_values_extractors[i]]
              .get_function_values(previous_solutions[p],
                                   previous_reactive_species_values[i][p]);
          }
      }

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_reactive_species.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            for (unsigned int i = 0; i < 4; i++)
              {
                // Shape functions for the phase order
                this->phi[i][q][k] =
                  this->fe_values_reactive_species[fe_values_extractors[i]]
                    .value(k, q);
                this->grad_phi[i][q][k] =
                  this->fe_values_reactive_species[fe_values_extractors[i]]
                    .gradient(k, q);
                this->hess_phi[i][q][k] =
                  this->fe_values_reactive_species[fe_values_extractors[i]]
                    .hessian(k, q);
                this->laplacian_phi[i][q][k] = trace(this->hess_phi[i][q][k]);
              }
          }
      }
  }


  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType           &current_solution,
                  const Parameters::ALE<dim> &ale)
  {
    FEValuesExtractors::Vector velocities(0);
    this->fe_values_fd.reinit(cell);

    this->fe_values_fd[velocities].get_function_values(current_solution,
                                                       velocity_values);

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
      }
  }

  /**
   * @brief Calculates the physical properties that may be required by the
   * Reactive species problem. Namely the surface tension and the Reactive
   * species mobility.
   */
  void
  calculate_physical_properties();

  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  std::vector<double>                  tracer_diffusivity;
  std::vector<double>                  tracer_diffusivity_0;
  std::vector<double>                  tracer_diffusivity_1;

  std::vector<FEValuesExtractors::Scalar> fe_values_extractors;

  // FEValues for the Reactive species problem
  FEValues<dim> fe_values_reactive_species;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Reactive species values
  std::vector<std::vector<double>>         reactive_species_values;
  std::vector<std::vector<Tensor<1, dim>>> reactive_species_gradients;
  std::vector<std::vector<double>>         reactive_species_laplacians;
  std::vector<std::vector<std::vector<double>>>
    previous_reactive_species_values;
  // TODO Initialize these vectors

  // Solid signed distance function
  std::vector<double> sdf_values;

  // Source term
  std::vector<std::vector<double>> source;
  // TODO Initialize this vectors

  // Shape functions
  std::vector<std::vector<std::vector<double>>>         phi;
  std::vector<std::vector<std::vector<Tensor<2, dim>>>> hess_phi;
  std::vector<std::vector<std::vector<double>>>         laplacian_phi;
  std::vector<std::vector<std::vector<Tensor<1, dim>>>> grad_phi;
  // TODO Initialize these vectors


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
};

#endif
