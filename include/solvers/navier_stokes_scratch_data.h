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


#ifndef lethe_navier_stokes_scratch_data_h
#define lethe_navier_stokes_scratch_data_h

#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/density_model.h>
#include <core/parameters.h>
#include <core/physical_property_model.h>
#include <core/rheological_model.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_filter.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/vof_filter.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;


/**
 * @brief Class that stores the information required by the assembly procedure
 * for a Navier-Stokes equation. Consequently, this class calculates
 * the velocity (values, gradients, laplacians) and the shape function
 * (values, gradients, laplacians) at all the gauss points for all degrees
 * of freedom and stores it into arrays. Additionally, the use can request
 * that this class gathers additional fields for physics which are coupled
 * to the Navier-Stokes equation, such as the VOF. This class
 * serves as a separation between the evaluation at the gauss point of the
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
   * evaluation, since this needs to be done at the cell level.
   *
   * @param simulation_control The SimulationControl object that holds
   * information related to the control of the steady-state or transient
   * simulation. This is used to extrapolate auxiliary physics solutions in time
   * for transient simulation.
   *
   * @param properties_manager The PhysicalPropertiesManager object that stores
   * physical property models.
   *
   * @param fe The FESystem used to solve the Navier-Stokes equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes
   * equations are solved
   *
   */
  NavierStokesScratchData(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const PhysicalPropertiesManager          &properties_manager,
    const FESystem<dim>                      &fe,
    const Quadrature<dim>                    &quadrature,
    const Mapping<dim>                       &mapping,
    const Quadrature<dim - 1>                &face_quadrature)
    : simulation_control(simulation_control)
    , properties_manager(properties_manager)
    , fe_values(mapping,
                fe,
                quadrature,
                update_values | update_quadrature_points | update_JxW_values |
                  update_gradients | update_hessians)
    , fe_face_values(mapping,
                     fe,
                     face_quadrature,
                     update_values | update_quadrature_points |
                       update_JxW_values | update_gradients | update_hessians |
                       update_normal_vectors)
  {
    allocate();

    // By default, the assembly of variables belonging to auxiliary physics is
    // disabled.
    gather_vof                               = false;
    gather_projected_phase_fraction_gradient = false;
    gather_curvature                         = false;
    gather_void_fraction                     = false;
    gather_particles_information             = false;
    gather_temperature                       = false;
    gather_cahn_hilliard                     = false;
    gather_hessian = properties_manager.is_non_newtonian();
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd The scratch data to be copied
   */
  NavierStokesScratchData(const NavierStokesScratchData<dim> &sd)
    : simulation_control(sd.simulation_control)
    , properties_manager(sd.properties_manager)
    , fe_values(sd.fe_values.get_mapping(),
                sd.fe_values.get_fe(),
                sd.fe_values.get_quadrature(),
                update_values | update_quadrature_points | update_JxW_values |
                  update_gradients | update_hessians)
    , fe_face_values(sd.fe_face_values.get_mapping(),
                     sd.fe_face_values.get_fe(),
                     sd.fe_face_values.get_quadrature(),
                     update_values | update_quadrature_points |
                       update_JxW_values | update_gradients | update_hessians |
                       update_normal_vectors)
  {
    allocate();

    // By default, the assembly of variables belonging to auxiliary physics is
    // disabled.
    gather_vof                               = false;
    gather_projected_phase_fraction_gradient = false;
    gather_curvature                         = false;
    gather_void_fraction                     = false;
    gather_particles_information             = false;
    gather_temperature                       = false;
    gather_cahn_hilliard                     = false;
    gather_hessian = properties_manager.is_non_newtonian();


    if (sd.gather_vof)
      enable_vof(sd.fe_values_vof->get_fe(),
                 sd.fe_values_vof->get_quadrature(),
                 sd.fe_values_vof->get_mapping(),
                 sd.filter);
    if (sd.gather_projected_phase_fraction_gradient)
      enable_projected_phase_fraction_gradient(
        sd.fe_values_projected_phase_fraction_gradient->get_fe(),
        sd.fe_values_projected_phase_fraction_gradient->get_quadrature(),
        sd.fe_values_projected_phase_fraction_gradient->get_mapping());
    if (sd.gather_curvature)
      enable_curvature(sd.fe_values_curvature->get_fe(),
                       sd.fe_values_curvature->get_quadrature(),
                       sd.fe_values_curvature->get_mapping());

    if (sd.gather_void_fraction)
      enable_void_fraction(sd.fe_values_void_fraction->get_fe(),
                           sd.fe_values_void_fraction->get_quadrature(),
                           sd.fe_values_void_fraction->get_mapping());
    if (sd.gather_particles_information)
      enable_particle_fluid_interactions(sd.max_number_of_particles_per_cell,
                                         sd.interpolated_void_fraction);
    if (sd.gather_temperature)
      enable_heat_transfer(sd.fe_values_temperature->get_fe(),
                           sd.fe_values_temperature->get_quadrature(),
                           sd.fe_values_temperature->get_mapping());
    if (sd.gather_cahn_hilliard)
      enable_cahn_hilliard(sd.fe_values_cahn_hilliard->get_fe(),
                           sd.fe_values_cahn_hilliard->get_quadrature(),
                           sd.fe_values_cahn_hilliard->get_mapping(),
                           sd.cahn_hilliard_filter);

    gather_hessian = sd.gather_hessian;
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
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for [u,p]
   *
   * @param previous_solutions The solutions at the previous time steps
   *
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
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         std::shared_ptr<Function<dim>> forcing_function,
         Tensor<1, dim>                 beta_force,
         const double                   pressure_scaling_factor)
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

    for (const unsigned int k : fe_values.dof_indices())
      {
        components[k] = fe.system_to_component_index(k).first;
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
    if (gather_hessian)
      this->fe_values[velocities].get_function_hessians(
        current_solution, this->velocity_hessians);

    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->velocity_divergences[q] = trace(this->velocity_gradients[q]);
      }

    // Gather pressure (values, gradient)
    fe_values[pressure].get_function_values(current_solution,
                                            this->pressure_values);
    fe_values[pressure].get_function_gradients(current_solution,
                                               this->pressure_gradients);
    this->pressure_scaling_factor = pressure_scaling_factor;

    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values[velocities].get_function_values(
          previous_solutions[p], previous_velocity_values[p]);
      }

    // Only gather the pressure when a pressure history is necessary
    // (compressible Navier-Stokes)
    if (!this->properties_manager.density_is_constant())
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        {
          this->fe_values[pressure].get_function_values(
            previous_solutions[p], previous_pressure_values[p]);
        }

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

    is_boundary_cell = cell->at_boundary();
    if (is_boundary_cell)
      {
        n_faces          = cell->n_faces();
        is_boundary_face = std::vector<bool>(n_faces, false);
        n_faces_q_points = fe_face_values.get_quadrature().size();
        boundary_face_id = std::vector<unsigned int>(n_faces);

        face_JxW = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));


        // Velocity and pressure values
        // First vector is face number, second quadrature point
        this->face_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        this->face_velocity_divergences = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        this->face_velocity_gradients =
          std::vector<std::vector<Tensor<2, dim>>>(
            n_faces, std::vector<Tensor<2, dim>>(n_faces_q_points));

        this->face_velocity_laplacians =
          std::vector<std::vector<Tensor<1, dim>>>(
            n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        this->face_pressure_values = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        this->face_pressure_gradients =
          std::vector<std::vector<Tensor<1, dim>>>(
            n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        this->face_normal = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        this->face_quadrature_points = std::vector<std::vector<Point<dim>>>(
          n_faces, std::vector<Point<dim>>(n_faces_q_points));

        this->face_div_phi_u = std::vector<std::vector<std::vector<double>>>(
          n_faces,
          std::vector<std::vector<double>>(n_faces_q_points,
                                           std::vector<double>(n_dofs)));

        this->face_phi_u =
          std::vector<std::vector<std::vector<Tensor<1, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<1, dim>>>(
              n_faces_q_points, std::vector<Tensor<1, dim>>(n_dofs)));

        this->face_hess_phi_u =
          std::vector<std::vector<std::vector<Tensor<3, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<3, dim>>>(
              n_faces_q_points, std::vector<Tensor<3, dim>>(n_dofs)));

        this->face_laplacian_phi_u =
          std::vector<std::vector<std::vector<Tensor<1, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<1, dim>>>(
              n_faces_q_points, std::vector<Tensor<1, dim>>(n_dofs)));

        this->face_grad_phi_u =
          std::vector<std::vector<std::vector<Tensor<2, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<2, dim>>>(
              n_faces_q_points, std::vector<Tensor<2, dim>>(n_dofs)));

        this->face_phi_p = std::vector<std::vector<std::vector<double>>>(
          n_faces,
          std::vector<std::vector<double>>(n_faces_q_points,
                                           std::vector<double>(n_dofs)));

        this->face_grad_phi_p =
          std::vector<std::vector<std::vector<Tensor<1, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<1, dim>>>(
              n_faces_q_points, std::vector<Tensor<1, dim>>(n_dofs)));
        for (const auto face : cell->face_indices())
          {
            is_boundary_face[face] = cell->face(face)->at_boundary();
            if (is_boundary_face[face])
              {
                fe_face_values.reinit(cell, face);
                n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
                boundary_face_id[face] = cell->face(face)->boundary_id();
                // Shape functions
                // First vector is face number, second quadrature point, third
                // DOF

                // Gather velocity (values, gradient and laplacian)
                this->fe_face_values[velocities].get_function_values(
                  current_solution, this->face_velocity_values[face]);
                this->fe_face_values[velocities].get_function_gradients(
                  current_solution, this->face_velocity_gradients[face]);
                this->fe_face_values[velocities].get_function_laplacians(
                  current_solution, this->face_velocity_laplacians[face]);

                // Gather pressure (values, gradient)
                this->fe_values[pressure].get_function_values(
                  current_solution, this->face_pressure_values[face]);
                this->fe_values[pressure].get_function_gradients(
                  current_solution, this->face_pressure_gradients[face]);

                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  {
                    this->face_JxW[face][q] = this->fe_face_values.JxW(q);
                    this->face_normal[face][q] =
                      this->fe_face_values.normal_vector(q);
                    this->face_quadrature_points[face][q] =
                      this->fe_face_values.quadrature_point(q);
                    for (const unsigned int k : fe_face_values.dof_indices())
                      {
                        // Velocity
                        this->face_phi_u[face][q][k] =
                          this->fe_face_values[velocities].value(k, q);
                        this->face_div_phi_u[face][q][k] =
                          this->fe_face_values[velocities].divergence(k, q);
                        this->face_grad_phi_u[face][q][k] =
                          this->fe_face_values[velocities].gradient(k, q);
                        this->face_hess_phi_u[face][q][k] =
                          this->fe_face_values[velocities].hessian(k, q);
                        for (int d = 0; d < dim; ++d)
                          this->face_laplacian_phi_u[face][q][k][d] =
                            trace(this->face_hess_phi_u[face][q][k][d]);
                        // Pressure
                        this->face_phi_p[face][q][k] =
                          this->fe_face_values[pressure].value(k, q);
                        this->face_grad_phi_p[face][q][k] =
                          this->fe_face_values[pressure].gradient(k, q);
                      }
                  }
              }
          }
      }
  }


  /**
   * @brief enable_vof Enables the collection of the VOF data by the scratch
   *
   * @param fe FiniteElement associated with the VOF.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param phase_filter_parameters Parameters for phase fraction filtering
   */

  void
  enable_vof(const FiniteElement<dim>          &fe,
             const Quadrature<dim>             &quadrature,
             const Mapping<dim>                &mapping,
             const Parameters::VOF_PhaseFilter &phase_filter_parameters);

  /**
   * @brief enable_vof Enables the collection of the VOF data by the scratch - function overload used in the copy constructor of NavierStokesScratchData
   *
   * @param fe FiniteElement associated with the VOF.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param filter Filter that is applied on the phase fraction
   */

  void
  enable_vof(const FiniteElement<dim>                       &fe,
             const Quadrature<dim>                          &quadrature,
             const Mapping<dim>                             &mapping,
             const std::shared_ptr<VolumeOfFluidFilterBase> &filter);

  void
  enable_projected_phase_fraction_gradient(
    const FiniteElement<dim> &fe_projected_phase_fraction_gradient,
    const Quadrature<dim>    &quadrature,
    const Mapping<dim>       &mapping);

  void
  enable_curvature(const FiniteElement<dim> &fe_curvature,
                   const Quadrature<dim>    &quadrature,
                   const Mapping<dim>       &mapping);

  /** @brief Reinitialize the content of the scratch for the vof
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the VOF FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [alpha]
   *
   * @param current_filtered_solution The present value of the solution for [alpha]_filtered
   *
   * @param previous_solutions The solutions at the previous time steps for [alpha]
   *
   *
   */

  template <typename VectorType>
  void
  reinit_vof(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType              &current_solution,
             const VectorType              &current_filtered_solution,
             const std::vector<VectorType> &previous_solutions)
  {
    this->fe_values_vof->reinit(cell);
    // Gather phase fraction (values, gradient)
    this->fe_values_vof->get_function_values(current_solution,
                                             this->phase_values);
    this->fe_values_vof->get_function_values(current_filtered_solution,
                                             this->filtered_phase_values);
    this->fe_values_vof->get_function_gradients(
      current_filtered_solution, this->filtered_phase_gradient_values);
    this->fe_values_vof->get_function_gradients(current_solution,
                                                this->phase_gradient_values);

    // Gather previous phase fraction values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_vof->get_function_values(previous_solutions[p],
                                                 previous_phase_values[p]);
      }
  }

  template <typename VectorType>
  void
  reinit_projected_phase_fraction_gradient(
    const typename DoFHandler<dim>::active_cell_iterator
                     &projected_phase_fraction_gradient_cell,
    const VectorType &current_projected_phase_fraction_gradient_solution)
  {
    this->fe_values_projected_phase_fraction_gradient->reinit(
      projected_phase_fraction_gradient_cell);

    FEValuesExtractors::Vector pfg(0);
    // Gather phase fraction gradient
    (*fe_values_projected_phase_fraction_gradient)[pfg].get_function_values(
      current_projected_phase_fraction_gradient_solution,
      this->projected_phase_fraction_gradient_values);
  }

  template <typename VectorType>
  void
  reinit_curvature(
    const typename DoFHandler<dim>::active_cell_iterator &curvature_cell,
    const VectorType &current_curvature_solution)
  {
    this->fe_values_curvature->reinit(curvature_cell);

    // Gather phase fraction gradient
    this->fe_values_curvature->get_function_values(current_curvature_solution,
                                                   this->curvature_values);
  }

  /**
   * @brief enable_void_fraction Enables the collection of the void fraction
   * data by the scratch
   *
   * @param fe FiniteElement associated with the void fraction
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */

  void
  enable_void_fraction(const FiniteElement<dim> &fe,
                       const Quadrature<dim>    &quadrature,
                       const Mapping<dim>       &mapping);

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
   */

  template <typename VectorType>
  void
  reinit_void_fraction(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const std::vector<VectorType>                        &previous_solutions)
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
    const unsigned int n_global_max_particles_per_cell,
    const bool         enable_void_fraction_interpolation);

  /** @brief Calculates the content of the scratch for the particle fluid
   * interactions
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the void fraction FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [epsilon]
   *
   * @param previous_solutions The solutions at the previous time steps for
   * [epsilon]
   *
   */

  template <typename VectorType>
  void
  reinit_particle_fluid_interactions(
    const typename DoFHandler<dim>::active_cell_iterator & /*cell*/,
    const VectorType & /*current_solution*/,
    const VectorType                      &previous_solution,
    const VectorType                      &void_fraction_solution,
    const Particles::ParticleHandler<dim> &particle_handler,
    DoFHandler<dim>                       &dof_handler,
    DoFHandler<dim>                       &void_fraction_dof_handler)
  {
    pic = particle_handler.particles_in_cell(this->fe_values.get_cell());

    average_particle_velocity = 0;
    cell_volume               = this->fe_values.get_cell()->measure();

    // Loop over particles in cell
    double total_particle_volume = 0;
    {
      unsigned int particle_i = 0;
      for (auto &particle : pic)
        {
          auto particle_properties = particle.get_properties();
          // Set the particle_fluid_interactions properties and vectors to 0
          for (int d = 0; d < dim; ++d)
            {
              particle_properties[DEM::PropertiesIndex::fem_force_x + d]  = 0.;
              particle_properties[DEM::PropertiesIndex::fem_torque_x + d] = 0.;
              undisturbed_flow_force[d]                                   = 0.;
            }

          // Stock the values of particle velocity in a tensor
          particle_velocity[particle_i][0] =
            particle_properties[DEM::PropertiesIndex::v_x];
          particle_velocity[particle_i][1] =
            particle_properties[DEM::PropertiesIndex::v_y];
          if constexpr (dim == 3)
            particle_velocity[particle_i][2] =
              particle_properties[DEM::PropertiesIndex::v_z];

          cell_void_fraction[particle_i] = 0;
          if (!interpolated_void_fraction)
            total_particle_volume +=
              M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
              (2 * dim);

          average_particle_velocity += particle_velocity[particle_i];
          particle_i++;
        }
      number_of_particles = particle_i;
    }



    if (!this->interpolated_void_fraction)
      {
        double cell_void_fraction_bulk = 0;
        cell_void_fraction_bulk =
          (cell_volume - total_particle_volume) / cell_volume;

        for (unsigned int j = 0; j < number_of_particles; ++j)
          cell_void_fraction[j] = cell_void_fraction_bulk;
      }

    if (number_of_particles == 0)
      return;


    // Calculate the average particle velocity within
    // the cell
    average_particle_velocity = average_particle_velocity / number_of_particles;

    // Create local vector that will be used to spawn an in-situ quadrature to
    // interpolate at the location of the particles
    std::vector<Point<dim>> particle_reference_location;
    std::vector<double>     particle_weights(number_of_particles, 1);

    // Loop over particles in cell and cache their reference location
    for (auto &particle : pic)
      {
        // Store particle positions and weights
        // Reference location of the particle
        particle_reference_location.push_back(
          particle.get_reference_location());
      }

    // Create a quadrature for the Navier-Stokes equations that is based on the
    // particle reference location
    Quadrature<dim> q_local(particle_reference_location, particle_weights);
    FEValues<dim>   fe_values_local_particles(this->fe_values.get_fe(),
                                            q_local,
                                            update_gradients | update_values |
                                              update_hessians);

    // Reallocate memory for the fields to be interpolated at the particle
    // location This has to be done for every cell because deal.II expects the
    // size of the vector to be strictly equal to the number of points in the
    // quadrature
    fluid_velocity_at_particle_location.resize(number_of_particles);
    fluid_velocity_laplacian_at_particle_location.resize(number_of_particles);
    fluid_velocity_curls_at_particle_location_2d.resize(number_of_particles);
    fluid_velocity_curls_at_particle_location_3d.resize(number_of_particles);
    fluid_pressure_gradients_at_particle_location.resize(number_of_particles);

    // Evaluate the relevant information at the
    // quadrature points to do the interpolation.
    const auto &velocity_cell =
      typename DoFHandler<dim>::cell_iterator(*this->fe_values.get_cell(),
                                              &dof_handler);

    fe_values_local_particles.reinit(velocity_cell);

    // Calculate all fluid properties at the particle location
    fe_values_local_particles[velocities].get_function_values(
      previous_solution, fluid_velocity_at_particle_location);

    fe_values_local_particles[velocities].get_function_laplacians(
      previous_solution, fluid_velocity_laplacian_at_particle_location);

    if constexpr (dim == 2)
      {
        fe_values_local_particles[velocities].get_function_curls(
          previous_solution, fluid_velocity_curls_at_particle_location_2d);
      }
    else if constexpr (dim == 3)
      {
        fe_values_local_particles[velocities].get_function_curls(
          previous_solution, fluid_velocity_curls_at_particle_location_3d);
      }

    fe_values_local_particles[pressure].get_function_gradients(
      previous_solution, fluid_pressure_gradients_at_particle_location);

    // Create a quadrature for the void fraction that is based on the particle
    // reference location
    if (this->interpolated_void_fraction)
      {
        FEValues<dim> fe_values_particles_void_fraction(
          this->fe_values_void_fraction->get_fe(), q_local, update_values);

        const auto &void_fraction_dh_cell =
          typename DoFHandler<dim>::cell_iterator(
            *this->fe_values_void_fraction->get_cell(),
            &void_fraction_dof_handler);

        fe_values_particles_void_fraction.reinit(void_fraction_dh_cell);

        fe_values_particles_void_fraction.get_function_values(
          void_fraction_solution, cell_void_fraction);
      }

    // Relative velocity and particle Reynolds
    unsigned int particle_i                  = 0;
    average_fluid_particle_relative_velocity = 0;

    // TODO, get the real viscosity at the particle localtion
    double kinematic_viscosity =
      properties_manager.get_kinematic_viscosity_scale();

    for (auto &particle : pic)
      {
        auto particle_properties = particle.get_properties();
        fluid_particle_relative_velocity_at_particle_location[particle_i] =
          fluid_velocity_at_particle_location[particle_i] -
          particle_velocity[particle_i];
        average_fluid_particle_relative_velocity +=
          fluid_particle_relative_velocity_at_particle_location[particle_i];

        Re_particle[particle_i] =
          1e-3 +
          cell_void_fraction[particle_i] *
            fluid_particle_relative_velocity_at_particle_location[particle_i]
              .norm() *
            particle_properties[DEM::PropertiesIndex::dp] /
            (kinematic_viscosity + DBL_MIN);
        particle_i++;
      }

    average_fluid_particle_relative_velocity =
      average_fluid_particle_relative_velocity / particle_i;
  }


  /**
   * @brief enable_heat_transfer Enables the collection of the heat transfer
   * data (Temperature field) by the scratch
   *
   * @param fe FiniteElement associated with the heat transfer.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */

  void
  enable_heat_transfer(const FiniteElement<dim> &fe,
                       const Quadrature<dim>    &quadrature,
                       const Mapping<dim>       &mapping);


  /** @brief Reinitialize the content of the scratch for the heat transfer
   * auxiliary physic
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the heat transfer FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for temperature
   *
   * @param previous_solutions Vector of \f$n\f$ @p VectorType containers of
   * previous temperature solutions. \f$n\f$ depends on the BDF scheme selected
   * for time-stepping.
   *
   */

  template <typename VectorType>
  void
  reinit_heat_transfer(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const std::vector<VectorType>                        &previous_solutions)
  {
    this->fe_values_temperature->reinit(cell);

    // Gather temperature
    this->fe_values_temperature->get_function_values(current_solution,
                                                     this->temperature_values);

    // Gather temperature gradient
    this->fe_values_temperature->get_function_gradients(
      current_solution, this->temperature_gradients);

    // Extrapolate temperature and temperature gradient to t+dt using the BDF
    // if the simulation is transient
    const auto method = this->simulation_control->get_assembly_method();
    if (is_bdf(method))
      {
        // Gather previous temperature and temperature gradient values
        for (unsigned int p = 0; p < previous_solutions.size(); ++p)
          {
            fe_values_temperature->get_function_values(
              previous_solutions[p], this->previous_temperature_values[p]);
            fe_values_temperature->get_function_gradients(
              previous_solutions[p], this->previous_temperature_gradients[p]);
          }
        // Extrapolate temperature and temperature gradient
        std::vector<double> time_vector =
          this->simulation_control->get_simulation_times();
        bdf_extrapolate(time_vector,
                        previous_temperature_values,
                        number_of_previous_solutions(method),
                        temperature_values);
        bdf_extrapolate(time_vector,
                        previous_temperature_gradients,
                        number_of_previous_solutions(method),
                        temperature_gradients);
      }
  }

  /**
   * @brief enable_cahn_hilliard Enables the collection of the CahnHilliard data
   * by the scratch
   *
   * @param fe FiniteElement associated with the CahnHilliard physics
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */
  void
  enable_cahn_hilliard(
    const FiniteElement<dim>       &fe,
    const Quadrature<dim>          &quadrature,
    const Mapping<dim>             &mapping,
    const Parameters::CahnHilliard &cahn_hilliard_parameters);

  /**
   * @brief Enables the collection of the CahnHilliard data by the scratch
   *
   * @param fe FiniteElement associated with the CahnHilliard physics
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   */
  void
  enable_cahn_hilliard(
    const FiniteElement<dim>                      &fe,
    const Quadrature<dim>                         &quadrature,
    const Mapping<dim>                            &mapping,
    const std::shared_ptr<CahnHilliardFilterBase> &cahn_hilliard_filter);


  /** @brief Reinitialize the content of the scratch for CH
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the CH FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [phi]
   *
   */
  template <typename VectorType>
  void
  reinit_cahn_hilliard(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const VectorType &current_filtered_solution)
  {
    this->fe_values_cahn_hilliard->reinit(cell);
    this->phase_order.component        = 0;
    this->chemical_potential.component = 1;

    // Gather phase fraction (values, gradients)
    this->fe_values_cahn_hilliard->operator[](phase_order)
      .get_function_values(current_solution,
                           this->phase_order_cahn_hilliard_values);
    this->fe_values_cahn_hilliard->operator[](chemical_potential)
      .get_function_values(current_solution,
                           this->chemical_potential_cahn_hilliard_values);
    this->fe_values_cahn_hilliard->operator[](phase_order)
      .get_function_gradients(current_solution,
                              this->phase_order_cahn_hilliard_gradients);

    // Gather filtered phase fraction (values, gradients)
    this->fe_values_cahn_hilliard->operator[](phase_order)
      .get_function_values(current_filtered_solution,
                           this->filtered_phase_order_cahn_hilliard_values);
  }


  /** @brief Calculates the physical properties. This function calculates the
   * physical properties that may be required by the fluid dynamics problem.
   * Namely the kinematic viscosity and, when required, the density.
   */
  void
  calculate_physical_properties();

  // For auxiliary physics solution extrapolation
  const std::shared_ptr<SimulationControl> simulation_control;

  // Physical properties
  const PhysicalPropertiesManager      properties_manager;
  std::map<field, std::vector<double>> fields;
  std::vector<double>                  density;
  double                               density_ref;
  double                               density_psi;
  std::vector<double>                  dynamic_viscosity;
  std::vector<double>                  kinematic_viscosity;
  double                               kinematic_viscosity_scale;
  /// Values of the kinematic viscosity used in the SUPG and PSPG
  /// stabilizations.
  std::vector<double> kinematic_viscosity_for_stabilization;
  /// Values of the dynamic viscosity used in the SUPG and PSPG stabilizations.
  std::vector<double>              dynamic_viscosity_for_stabilization;
  std::vector<double>              thermal_expansion;
  std::vector<double>              grad_kinematic_viscosity_shear_rate;
  std::vector<std::vector<double>> previous_density;

  // Pressure scaling factor to facilitate different scales between velocity and
  // pressure
  double pressure_scaling_factor;

  // For VOF and CH simulations. Present properties for fluid 0 and 1.
  std::vector<double> density_0;
  std::vector<double> density_1;
  double              density_ref_0;
  double              density_ref_1;
  double              density_psi_0;
  double              density_psi_1;
  std::vector<double> compressibility_multiplier;
  std::vector<double> dynamic_viscosity_0;
  std::vector<double> dynamic_viscosity_1;
  std::vector<double> kinematic_viscosity_0;
  std::vector<double> kinematic_viscosity_1;
  /// Scale of the kinematic viscosity for fluid 0.
  double kinematic_viscosity_scale_0;
  /// Scale of the kinematic viscosity for fluid 1.
  double kinematic_viscosity_scale_1;
  /// Values of the dynamic viscosity used in the SUPG and PSPG stabilizations
  /// for fluid 0.
  std::vector<double> dynamic_viscosity_for_stabilization_0;
  /// Values of the dynamic viscosity used in the SUPG and PSPG stabilizations
  /// for fluid 1.
  std::vector<double> dynamic_viscosity_for_stabilization_1;
  std::vector<double> thermal_expansion_0;
  std::vector<double> thermal_expansion_1;
  std::vector<double> surface_tension;
  std::vector<double> surface_tension_gradient;

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

  // Components index
  std::vector<unsigned int> components;

  // Velocity and pressure values
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<double>                      velocity_divergences;
  std::vector<Tensor<2, dim>>              velocity_gradients;
  std::vector<Tensor<1, dim>>              velocity_laplacians;
  std::vector<Tensor<3, dim>>              velocity_hessians;
  std::vector<double>                      shear_rate;
  std::vector<double>                      pressure_values;
  std::vector<Tensor<1, dim>>              pressure_gradients;
  std::vector<std::vector<double>>         previous_pressure_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;

  // Shape functions
  std::vector<std::vector<double>>         div_phi_u;
  std::vector<std::vector<Tensor<1, dim>>> phi_u;
  std::vector<std::vector<Tensor<3, dim>>> hess_phi_u;
  std::vector<std::vector<Tensor<1, dim>>> laplacian_phi_u;
  std::vector<std::vector<Tensor<2, dim>>> grad_phi_u;
  std::vector<std::vector<double>>         phi_p;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_p;

  /**
   * Scratch component for the VOF auxiliary physics
   */
  bool                             gather_vof;
  unsigned int                     n_dofs_vof;
  std::vector<double>              phase_values;
  std::vector<double>              filtered_phase_values;
  std::vector<std::vector<double>> previous_phase_values;
  std::vector<Tensor<1, dim>>      filtered_phase_gradient_values;
  std::vector<Tensor<1, dim>>      phase_gradient_values;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>>           fe_values_vof;
  std::shared_ptr<VolumeOfFluidFilterBase> filter; // Phase fraction filter

  bool                           gather_projected_phase_fraction_gradient;
  bool                           gather_curvature;
  std::shared_ptr<FEValues<dim>> fe_values_projected_phase_fraction_gradient;
  std::shared_ptr<FEValues<dim>> fe_values_curvature;
  std::vector<Tensor<1, dim>>    projected_phase_fraction_gradient_values;
  std::vector<double>            curvature_values;

  /**
   * Scratch component for the void fraction auxiliary physics
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
  bool                        interpolated_void_fraction;
  std::vector<Tensor<1, dim>> particle_velocity;
  Tensor<1, dim>              average_particle_velocity;
  std::vector<Tensor<1, dim>> fluid_velocity_at_particle_location;
  std::vector<Tensor<1, dim>>
                 fluid_particle_relative_velocity_at_particle_location;
  Tensor<1, dim> average_fluid_particle_relative_velocity;
  std::vector<Tensor<1, dim>> fluid_pressure_gradients_at_particle_location;
  std::vector<Tensor<1, dim>> fluid_velocity_laplacian_at_particle_location;
  std::vector<Tensor<1, 1>>   fluid_velocity_curls_at_particle_location_2d;
  std::vector<Tensor<1, 3>>   fluid_velocity_curls_at_particle_location_3d;
  std::vector<Point<dim>>     particle_reference_location;
  std::vector<double>         particle_weights;
  std::vector<double>         cell_void_fraction;
  std::vector<double>         Re_particle;
  unsigned int                max_number_of_particles_per_cell;
  unsigned int                number_of_particles;
  typename Particles::ParticleHandler<dim>::particle_iterator_range pic;
  double                                                            cell_volume;
  double                                                            beta_drag;
  Tensor<1, dim> undisturbed_flow_force;

  /**
   * Scratch component for the heat transfer
   */
  bool                                     gather_temperature;
  unsigned int                             n_dofs_heat_transfer;
  std::vector<double>                      temperature_values;
  std::vector<std::vector<double>>         previous_temperature_values;
  std::vector<Tensor<1, dim>>              temperature_gradients;
  std::vector<std::vector<Tensor<1, dim>>> previous_temperature_gradients;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_temperature;

  /**
   * Scratch component for the CahnHilliard auxiliary physics
   */
  bool                        gather_cahn_hilliard;
  unsigned int                n_dofs_cahn_hilliard;
  std::vector<double>         phase_order_cahn_hilliard_values;
  std::vector<double>         filtered_phase_order_cahn_hilliard_values;
  std::vector<Tensor<1, dim>> phase_order_cahn_hilliard_gradients;
  std::vector<double>         chemical_potential_cahn_hilliard_values;

  std::shared_ptr<CahnHilliardFilterBase>
    cahn_hilliard_filter; // Phase order fraction filter

  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_cahn_hilliard;
  FEValuesExtractors::Scalar     phase_order;
  FEValuesExtractors::Scalar     chemical_potential;

  /**
   * Is boundary cell indicator
   */
  bool is_boundary_cell;

  // If a rheological model is being used for a non-Newtonian flow
  bool gather_hessian;

  FEFaceValues<dim> fe_face_values;

  unsigned int n_faces;
  unsigned int n_faces_q_points;
  unsigned int face_n_dofs;

  // If boundary cell indicator
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;

  // Quadrature
  std::vector<std::vector<double>>         face_JxW;
  std::vector<std::vector<Point<dim>>>     face_quadrature_points;
  std::vector<std::vector<Tensor<1, dim>>> face_normal;

  // Velocity and pressure values
  // First vector is face number, second quadrature point
  std::vector<std::vector<Tensor<1, dim>>> face_velocity_values;
  std::vector<std::vector<double>>         face_velocity_divergences;
  std::vector<std::vector<Tensor<2, dim>>> face_velocity_gradients;
  std::vector<std::vector<Tensor<1, dim>>> face_velocity_laplacians;
  std::vector<std::vector<double>>         face_pressure_values;
  std::vector<std::vector<Tensor<1, dim>>> face_pressure_gradients;

  // Shape functions
  // First vector is face number, second quadrature point, third DOF
  std::vector<std::vector<std::vector<double>>>         face_div_phi_u;
  std::vector<std::vector<std::vector<Tensor<1, dim>>>> face_phi_u;
  std::vector<std::vector<std::vector<Tensor<3, dim>>>> face_hess_phi_u;
  std::vector<std::vector<std::vector<Tensor<1, dim>>>> face_laplacian_phi_u;
  std::vector<std::vector<std::vector<Tensor<2, dim>>>> face_grad_phi_u;
  std::vector<std::vector<std::vector<double>>>         face_phi_p;
  std::vector<std::vector<std::vector<Tensor<1, dim>>>> face_grad_phi_p;
};

#endif
