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

#include <dem/dem.h>
#include <dem/dem_properties.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/particle_handler.h>


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
    gather_free_surface          = false;
    gather_void_fraction         = false;
    gather_particles_information = false;
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

    if (sd.gather_void_fraction)
      enable_void_fraction(sd.fe_values_void_fraction->get_fe(),
                           sd.fe_values_void_fraction->get_quadrature(),
                           sd.fe_values_void_fraction->get_mapping());
    if (sd.gather_particles_information)
      enable_particle_fluid_interactions(sd.max_number_of_particles_per_cell);
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

        const unsigned int component_mass =
          fe.system_to_component_index(dim).first;
        this->mass_source[q] = this->rhs_force[q](component_mass);

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
                      const Mapping<dim> &      mapping);

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

  /**
   * @brief enable_void_fraction Enables the collection of the void fraction data by the scratch
   *
   * @param fe FiniteElement associated with the void fraction
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */

  void
  enable_void_fraction(const FiniteElement<dim> &fe,
                       const Quadrature<dim> &   quadrature,
                       const Mapping<dim> &      mapping);

  /** @brief Reinitialize the content of the scratch for the void fraction
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the void fraction FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [epsilon]
   *
   * @param previous_solutions The solutions at the previous time steps for [epsilon]
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods) for [epsilon]
   *
   */

  template <typename VectorType>
  void
  reinit_void_fraction(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType &                                    current_solution,
    const std::vector<VectorType> &                       previous_solutions,
    const std::vector<VectorType> & /*solution_stages*/)
  {
    this->fe_values_void_fraction->reinit(cell);

    // Gather void fraction (values, gradient)
    this->fe_values_void_fraction->get_function_values(
      current_solution, this->void_fraction_values);
    this->fe_values_void_fraction->get_function_gradients(
      current_solution, this->void_fraction_gradient_values);

    // Gather previous void fraction fraction values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_void_fraction->get_function_values(
          previous_solutions[p], previous_void_fraction_values[p]);
      }
  }

  /**
   * @brief enable_particle_fluid_interactions Enables the calculation of the drag force by the scratch
   *
   * @param fe FiniteElement associated with the void fraction
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */

  void
  enable_particle_fluid_interactions(
    const unsigned int n_global_max_particles_per_cell);

  /** @brief Calculate the content of the scratch for the particle fluid interactions
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the void fraction FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [epsilon]
   *
   * @param previous_solutions The solutions at the previous time steps for [epsilon]
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods) for [epsilon]
   *
   */

  template <typename VectorType>
  void
  reinit_particle_fluid_interactions(
    const VectorType                       velocity_solution,
    const VectorType                       void_fraction_solution,
    const Particles::ParticleHandler<dim> &particle_handler,
    DoFHandler<dim> &                      dof_handler,
    DoFHandler<dim> &                      void_fraction_dof_handler)
  {
    const FiniteElement<dim> &fe = this->fe_values.get_fe();
    const FiniteElement<dim> &fe_void_fraction =
      this->fe_values_void_fraction->get_fe();

    const unsigned int                   dofs_per_cell = fe.dofs_per_cell;
    std::vector<types::global_dof_index> fluid_dof_indices(dofs_per_cell);
    std::vector<types::global_dof_index> void_fraction_dof_indices(
      fe_void_fraction.dofs_per_cell);

    pic = particle_handler.particles_in_cell(this->fe_values.get_cell());

    unsigned int particle_number = 0;
    average_particle_velocity    = 0;
    cell_volume                  = this->fe_values.get_cell()->measure();

    // Get the local dof indices for velocity and void fraction field
    // interpolation at particle's position
    const auto &dh_cell =
      typename DoFHandler<dim>::cell_iterator(*this->fe_values.get_cell(),
                                              &dof_handler);
    dh_cell->get_dof_indices(fluid_dof_indices);

    const auto &void_fraction_dh_cell = typename DoFHandler<dim>::cell_iterator(
      *this->fe_values_void_fraction->get_cell(), &void_fraction_dof_handler);
    void_fraction_dh_cell->get_dof_indices(void_fraction_dof_indices);

    // Loop over particles in cell
    for (auto &particle : pic)
      {
        auto particle_properties = particle.get_properties();

        cell_void_fraction[particle_number]                  = 0;
        fluid_velocity_at_particle_location[particle_number] = 0;

        // Stock the values of particle velocity in a tensor
        particle_velocity[particle_number][0] =
          particle_properties[DEM::PropertiesIndex::v_x];
        particle_velocity[particle_number][1] =
          particle_properties[DEM::PropertiesIndex::v_y];
        if (dim == 3)
          particle_velocity[particle_number][2] =
            particle_properties[DEM::PropertiesIndex::v_z];

        // Interpolate velocity and void fraction at particle position
        // Reference location of the particle
        auto reference_location = particle.get_reference_location();

        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            const auto comp_j = fe.system_to_component_index(j).first;
            if (comp_j < dim)
              {
                auto &evaluation_point = velocity_solution;
                fluid_velocity_at_particle_location[particle_number][comp_j] +=
                  evaluation_point[fluid_dof_indices[j]] *
                  fe.shape_value(j, reference_location);
              }
          }
        for (unsigned int j = 0; j < fe_void_fraction.dofs_per_cell; ++j)
          {
            cell_void_fraction[particle_number] +=
              void_fraction_solution[void_fraction_dof_indices[j]] *
              fe_void_fraction.shape_value(j, reference_location);
          }

        average_particle_velocity += particle_velocity[particle_number];
        particle_number += 1;
      }

    if (particle_number != 0)
      {
        average_particle_velocity = average_particle_velocity / particle_number;
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
  std::vector<double>         mass_source;

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


  /**
   * Scratch component for the void fractoin auxiliary physics
   */
  bool                             gather_void_fraction;
  unsigned int                     n_dofs_void_fraction;
  std::vector<double>              void_fraction_values;
  std::vector<std::vector<double>> previous_void_fraction_values;
  std::vector<Tensor<1, dim>>      void_fraction_gradient_values;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_void_fraction;

  /**
   * Scratch component for the particle fluid interaction auxiliary physics
   */
  bool                        gather_particles_information;
  std::vector<Tensor<1, dim>> particle_velocity;
  Tensor<1, dim>              average_particle_velocity;
  std::vector<Tensor<1, dim>> fluid_velocity_at_particle_location;
  std::vector<double>         cell_void_fraction;
  unsigned int                max_number_of_particles_per_cell;
  unsigned int                number_of_locally_owned_particles;
  typename Particles::ParticleHandler<dim>::particle_iterator_range pic;
  double                                                            cell_volume;
  double                                                            beta_drag;
  std::vector<Tensor<1, dim>> fluid_particle_force;
  std::vector<unsigned int>   local_particle_id;
  unsigned int                particle_index;
};

#endif
