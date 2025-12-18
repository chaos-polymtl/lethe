// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_navier_stokes_scratch_data_h
#define lethe_navier_stokes_scratch_data_h

#include "core/parameters_cfd_dem.h"
#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/parameters.h>
#include <core/physical_property_model.h>
#include <core/sdirk_stage_data.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_filter.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>
#include <solvers/vof_filter.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_system.h>

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
 *
 * @ingroup solvers
 **/

template <int dim>
class NavierStokesScratchData : public PhysicsScratchDataBase
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
    , sdirk_table(
        time_stepping_is_sdirk(simulation_control->get_assembly_method()) ?
          ::sdirk_table(simulation_control->get_assembly_method()) :
          SDIRKTable())
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
    gather_mortar                            = false;
    gather_particle_field_project            = false;
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
    , sdirk_table(sd.sdirk_table)
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
    gather_mortar                            = false;
    gather_particle_field_project            = false;
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
    if (sd.gather_mortar)
      enable_mortar();

    if (sd.gather_particle_field_project)
      enable_particle_field_projection(
        sd.fe_values.get_quadrature(),
        sd.fe_values.get_mapping(),
        sd.fe_values_particle_drag->get_fe(),
        sd.fe_values_particle_two_way_coupling_force->get_fe(),
        sd.fe_values_particle_velocity->get_fe(),
        sd.fe_values_particle_momentum_transfer_coefficient->get_fe());

    gather_hessian = sd.gather_hessian;
  }

  /**
   * @brief Allocates the memory for the scratch
   *
   * This function allocates the necessary memory for all members of the scratch
   */
  void
  allocate() override;

  /**
   * @brief Reinitializes the content of the scratch.
   *
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the FE which is used to fill the
   * FeValues.
   *
   * @param[in] current_solution The present value of the solution for [u,p].
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   * @param[in] forcing_function The function describing the momentum/mass
   * source term.
   *
   * @param[in] beta_force The additional force for flow control.
   * TODO : Deprecate this argument and pass it to the constructor of the
   * assembler
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         const VectorType              &sum_over_previous_stages,
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

    // Compute cell diameter
    double cell_measure =
      compute_cell_measure_with_JxW(this->fe_values.get_JxW_values());
    this->cell_size = compute_cell_diameter<dim>(cell_measure, fe.degree);

    // For the SDIRK methods, \sum_{j=1}^{i-1} a_{ij} k_j is needed for the
    // assembler
    if (this->simulation_control->is_sdirk())
      {
        this->fe_values[velocities].get_function_values(
          sum_over_previous_stages, this->sdirk_stage_sum);
      }

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

    // Gather velocity for stabilization (same as velocity_values unless ALE is
    // enabled. The ALE correction is made within the correponding reinit
    // function)
    this->velocity_for_stabilization = this->velocity_values;

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

        face_JxW.reinit(n_faces, n_faces_q_points);

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

        this->face_normal.reinit(n_faces, n_faces_q_points);

        this->face_quadrature_points = std::vector<std::vector<Point<dim>>>(
          n_faces, std::vector<Point<dim>>(n_faces_q_points));

        this->face_div_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

        this->face_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);
        this->face_hess_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

        this->face_laplacian_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

        this->face_grad_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);
        this->face_phi_p.reinit(n_faces, n_faces_q_points, n_dofs);
        this->face_grad_phi_p.reinit(n_faces, n_faces_q_points, n_dofs);

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
                this->fe_face_values[pressure].get_function_values(
                  current_solution, this->face_pressure_values[face]);
                this->fe_face_values[pressure].get_function_gradients(
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

  /**
   * @brief Reinitialize the content of the scratch for the vof
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
   */

  template <typename VectorType>
  void
  reinit_vof(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType              &current_solution,
             const VectorType              &current_filtered_solution,
             const std::vector<VectorType> &previous_solutions)
  {
    Assert(
      gather_vof,
      ExcMessage(
        "You are trying to reinit the VOF model in a cell, but you did not enable VOF for the scratch data (gather_vof=false)."));

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
    Assert(
      gather_curvature,
      ExcMessage(
        "You are trying to reinit the curvature in a cell, but you did not enable curvature for the scratch data (gather_curvature=false)."));

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

  /**
   * @brief enable_particle_field_projection Enables the collection of the particle
   * fields projection data by the scratch
   *
   * @param[in] quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param[in] mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param[in] fe_particle_drag_proj FiniteElement associated with the
   * projected particle drag force
   *
   * @param[in] fe_particle_two_way_coupling_force_proj FiniteElement associated
   * with the projected particle two-way coupling force
   *
   * @param[in] fe_particle_velocity_proj FiniteElement associated with the
   * projected particle velocity
   *
   * @param[in] fe_particle_momentum_transfer_coefficient FiniteElement
   * associated with the projected particle momentum transfer coefficient
   */
  void
  enable_particle_field_projection(
    const Quadrature<dim>    &quadrature,
    const Mapping<dim>       &mapping,
    const FiniteElement<dim> &fe_particle_drag_proj,
    const FiniteElement<dim> &fe_particle_two_way_coupling_force_proj,
    const FiniteElement<dim> &fe_particle_velocity_proj,
    const FiniteElement<dim> &fe_particle_momentum_transfer_coefficient);

  /**
   *  @brief Reinitialize the content of the scratch for the void fraction
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the void fraction FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [epsilon]
   *
   * @param previous_solutions The solutions at the previous time steps for [epsilon]
   */

  template <typename VectorType>
  void
  reinit_void_fraction(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const std::vector<VectorType>                        &previous_solutions)
  {
    Assert(
      gather_void_fraction,
      ExcMessage(
        "You are trying to reinit the void fraction in a cell, but you did not enable void fraction for the scratch data (gather_void_fraction=false)."));

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
   * @brief enable_particle_fluid_interactions Enables the calculation of the fluid information at the particle location for the scratch.
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

  /**
   * @brief Reinitializes the fluid forces and torques on the particles in the cell to zero
   */
  void
  reinit_particle_fluid_forces();

  /**
   * @brief Extracts the velocity of the particles and calculates their total volume
   * in the cell. Both outcomes are stored as members of the scratch data class.
   */
  void
  extract_particle_properties();

  /**
   * @brief Computes the cell void fraction
   *
   * @param[in] total_particle_volume The total volume of the particles in the
   * cell
   */

  void
  calculate_cell_void_fraction(const double &total_particle_volume);

  /**
   * @brief Creates an object of type Quadrature<dim> that contains the
   * particle reference location. This object can be used to interpolate any
   * function known at the dofs at the location of the particles.
   *
   * @return Quadrature type object that contains the location of the particles
   * relative to the cell's frame of reference.
   */

  Quadrature<dim>
  gather_particles_reference_location();


  /**
   * @brief Interpolates the velocity and pressure of the fluid, as well as the
   * pressure gradient, and the laplacian, curl and gradient of the velocity,
   * at the locations of the particles.
   *
   * @param[in] q_particles_location Quadrature type object that contains the
   * location of the particles relative to the cell's frame of reference.
   *
   * @param[in] velocity_cell The active cell associated with the velocity and
   * pressure DoFHandler
   *
   * @param[in] present_velocity_pressure_solution The solution (velocity and
   * pressure) at the current time step. This solution is used in the implicit
   * coupling to establish the particle-fluid force.
   *
   * @param[in] previous_velocity_pressure_solution The solution (velocity and
   * pressure) at the previous time step. This solution is used in the explicit
   * and semi-implicit coupling to calculate the particle-fluid force.
   *
   * @param[in] drag_coupling Indicator for the type of coupling that is used.
   * This parameter essentially decides which velocity_pressure solution is
   * interpolated at the location of the particles.
   */

  template <typename VectorType>
  void
  calculate_fluid_fields_at_particle_location(
    const Quadrature<dim>                                &q_particles_location,
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const VectorType               &present_velocity_pressure_solution,
    const VectorType               &previous_velocity_pressure_solution,
    const Parameters::DragCoupling &drag_coupling)
  {
    FEValues<dim> fe_values_local_particles(this->fe_values.get_fe(),
                                            q_particles_location,
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

    // Take velocity_pressure_solution according to the type of coupling used.
    const auto &velocity_pressure_solution =
      drag_coupling == Parameters::DragCoupling::fully_implicit ?
        present_velocity_pressure_solution :
        previous_velocity_pressure_solution;

    fe_values_local_particles.reinit(velocity_cell);

    // Calculate all fluid properties at the particle location
    fe_values_local_particles[velocities].get_function_values(
      velocity_pressure_solution, fluid_velocity_at_particle_location);

    fe_values_local_particles[velocities].get_function_laplacians(
      velocity_pressure_solution,
      fluid_velocity_laplacian_at_particle_location);

    if constexpr (dim == 2)
      {
        fe_values_local_particles[velocities].get_function_curls(
          velocity_pressure_solution,
          fluid_velocity_curls_at_particle_location_2d);
      }
    else if constexpr (dim == 3)
      {
        fe_values_local_particles[velocities].get_function_curls(
          velocity_pressure_solution,
          fluid_velocity_curls_at_particle_location_3d);
      }

    fe_values_local_particles[pressure].get_function_gradients(
      velocity_pressure_solution,
      fluid_pressure_gradients_at_particle_location);
  }

  /**
   * @brief Interpolates the void fraction at the locations of the particles.
   *
   * @param[in] q_particles_location Quadrature type object that contains the
   * location of the particles relative to the cell's frame of reference.
   *
   * @param[in] void_fraction_cell The active cell associated with the void
   * fraction DoFHandler
   *
   * @param[in] void_fraction_solution The void fraction calculated with one of
   * the methods of the VoidFractionBase class.
   */

  template <typename VectorType>
  void
  calculate_void_fraction_at_particle_location(
    const Quadrature<dim>                                &q_particles_location,
    const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
    const VectorType &void_fraction_solution)
  {
    Assert(
      gather_void_fraction,
      ExcMessage(
        "gather_void_fraction has been set to false, yet you are trying to gather the void fraction at the location of the particles. The scratch data is currently unaware of the finite element interpolation for the void fraction and the simulation will abort."));

    FEValues<dim> fe_values_particles_void_fraction(
      this->fe_values_void_fraction->get_fe(),
      q_particles_location,
      update_values);

    Assert(
      cell_void_fraction.size() == q_particles_location.size(),
      ExcMessage(
        "The vector for the void fraction at the particle location does not have the same size as the quadrature used to evaluate it."));
    fe_values_particles_void_fraction.reinit(void_fraction_cell);

    fe_values_particles_void_fraction.get_function_values(
      void_fraction_solution, cell_void_fraction);
  }

  /**
   * @brief Calculates the properties of the fluid at the locations of the particles.
   * At the moment, only constant properties within the same fluid are
   * supported. When two fluids are present and VOF is used, the properties are
   * calculated based on the filtered VOF solution interpolated at the
   * location of the particles. These properties are used in the forces
   * calculations in the VANS equations.
   */

  void
  calculate_fluid_properties_at_particle_location()
  {
    AssertThrow(
      properties_manager.density_is_constant(),
      ExcMessage(
        "The lethe-fluid-particles solver only supports constant densities"));
    AssertThrow(
      !properties_manager.get_rheology()
          ->is_non_newtonian_rheological_model() &&
        !properties_manager.get_rheology()
           ->is_non_newtonian_rheological_model(),
      ExcMessage(
        "The lethe-fluid-particles solver only supports constant rheology model"));

    if (gather_vof)
      {
        for (unsigned int i_particle = 0; i_particle < number_of_particles;
             ++i_particle)
          {
            density_at_particle_location[i_particle] = calculate_point_property(
              filtered_phase_values_at_particle_location[i_particle],
              this->density_ref_0,
              this->density_ref_1);
            double dynamic_viscosity_at_particle_location =
              calculate_point_property(
                filtered_phase_values_at_particle_location[i_particle],
                this->kinematic_viscosity_scale_0 * this->density_ref_0,
                this->kinematic_viscosity_scale_1 * this->density_ref_1);

            kinematic_viscosity_at_particle_location[i_particle] =
              dynamic_viscosity_at_particle_location /
              density_at_particle_location[i_particle];
          }

        // Properties have been calculated, return
        return;
      }

    // Regular case without VOF
    for (unsigned int i_particle = 0; i_particle < number_of_particles;
         ++i_particle)
      {
        // Gather the kinematic viscosity and density at the particle
        // location assuming a constant kinematic viscosity and density in a
        // single fluid
        kinematic_viscosity_at_particle_location[i_particle] =
          kinematic_viscosity_scale;
        density_at_particle_location[i_particle] = density_scale;
      }
  }

  /**
   * @brief Calculates the velocity of the fluid relative to that of the particle and
   * the particle Reynolds number at the location of the particles. These will
   * be used in forces calculations in the vans equations.
   */

  void
  calculate_force_parameters_at_particle_location()
  {
    // Tolerance for the calculation of the particle Reynolds number.
    // TODO -> Revisit if that is not a too high value.
    const double Re_particle_tolerance       = 1e-3;
    unsigned int i_particle                  = 0;
    average_fluid_particle_relative_velocity = 0;

    for (auto &particle : pic)
      {
        auto particle_properties = particle.get_properties();

        fluid_particle_relative_velocity_at_particle_location[i_particle] =
          fluid_velocity_at_particle_location[i_particle] -
          particle_velocity[i_particle];
        average_fluid_particle_relative_velocity +=
          fluid_particle_relative_velocity_at_particle_location[i_particle];

        Re_particle[i_particle] =
          Re_particle_tolerance +
          cell_void_fraction[i_particle] *
            fluid_particle_relative_velocity_at_particle_location[i_particle]
              .norm() *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] /
            (kinematic_viscosity_at_particle_location[i_particle] + DBL_MIN);
        i_particle++;
      }

    average_fluid_particle_relative_velocity =
      average_fluid_particle_relative_velocity / i_particle;
  }

  /**
   * @brief Interpolates the filtered VOF solution at the location of the particles.
   * The latter values are used in calculating the density and viscosity of the
   * fluid at the particles' locations when VOF is used.
   *
   * @param[in] q_particles_location Quadrature type object that contains the
   * location of the particles relative to the cell's frame of reference.
   *
   * @param[in] phase_cell The active cell associated with the VOF DoFHandler
   *
   * @param[in] current_filtered_solution The present value of the filtered VOF
   * solution
   */

  template <typename VectorType>
  void
  calculate_vof_at_particle_location(
    const Quadrature<dim>                                &q_particles_location,
    const typename DoFHandler<dim>::active_cell_iterator &phase_cell,
    const VectorType &current_filtered_solution)
  {
    Assert(
      gather_vof,
      ExcMessage(
        "gather_vof has been set to false, yet you are trying to gather the VOF at the location of the particles. The scratch data is currently unaware of the finite element interpolation for VOF and the simulation will abort."));


    FEValues<dim> fe_values_vof_local_particles((*this->fe_values_vof).get_fe(),
                                                q_particles_location,
                                                update_values |
                                                  update_quadrature_points |
                                                  update_JxW_values);

    filtered_phase_values_at_particle_location.resize(number_of_particles);

    fe_values_vof_local_particles.reinit(phase_cell);

    fe_values_vof_local_particles.get_function_values(
      current_filtered_solution, filtered_phase_values_at_particle_location);
  }

  /**
   * @brief Calculates the variables needed to compute the particle fluid interactions
   * in the VANS equations.
   *
   * @param[in] velocity_cell The active cell associated with the velocity and
   * pressure DoFHandler
   *
   * @param[in] void_fraction_cell The active cell associated with the void
   * fraction DoFHandler
   *
   * @param[in] present_velocity_pressure_solution The solution (velocity and
   * pressure) at the current time step. This solution is used in the implicit
   * coupling to establish the particle-fluid force.
   *
   * @param[in] previous_velocity_pressure_solution The solution (velocity and
   * pressure) at the previous time step. This solution is used in the explicit
   * and semi-implicit coupling to calculate the particle-fluid force.
   *
   * @param[in] void_fraction_solution The void fraction value calculated with
   * one of the methods of the VoidFractionBase class.
   *
   * @param[in] particle_handler The particle handler object that stores and
   * manages the particles in the simulations
   *
   * @param[in] drag_coupling Indicator for the type of coupling that is used.
   * This parameter essentially decides which velocity_pressure solution is
   * interpolated at the location of the particles.
   */

  template <typename VectorType>
  void
  reinit_particle_fluid_interactions(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
    const VectorType                      &present_velocity_pressure_solution,
    const VectorType                      &previous_velocity_pressure_solution,
    const VectorType                      &void_fraction_solution,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::DragCoupling        &drag_coupling)
  {
    Assert(
      gather_particles_information,
      ExcMessage(
        "You are trying to reinit the fluid information at the particle location within a cell,"
        " but you did not enable the gathering of the information at the particle location for the scratch data (gather_particles_information=false)."));

    pic = particle_handler.particles_in_cell(velocity_cell);

    extract_particle_properties();
    reinit_particle_fluid_forces();

    calculate_cell_void_fraction(total_particle_volume);

    if (number_of_particles == 0)
      return;

    Quadrature<dim> q_particles_location =
      gather_particles_reference_location();

    calculate_fluid_fields_at_particle_location(
      q_particles_location,
      velocity_cell,
      present_velocity_pressure_solution,
      previous_velocity_pressure_solution,
      drag_coupling);

    if (this->interpolated_void_fraction)
      {
        calculate_void_fraction_at_particle_location(q_particles_location,
                                                     void_fraction_cell,
                                                     void_fraction_solution);
      }
    calculate_fluid_properties_at_particle_location();
    calculate_force_parameters_at_particle_location();
  }

  /**
   * @brief Calculates the variables needed to compute the particle fluid interactions
   * in the VANS equations.This version of the function is used when VOF is used
   * with CFD-DEM.
   *
   * @param[in] velocity_cell The active cell associated with the velocity and
   * pressure DoFHandler.
   *
   * @param[in] void_fraction_cell The active cell associated with the void
   * fraction DoFHandler.
   *
   * @param[in] phase_cell The active cell associated with the VOF DoFHandler.
   *
   * @param[in] present_velocity_pressure_solution The solution (velocity and
   * pressure) at the current time step. This solution is used in the implicit
   * coupling to establish the particle-fluid force.
   *
   * @param[in] previous_velocity_pressure_solution The solution (velocity and
   * pressure) at the previous time step. This solution is used in the explicit
   * and semi-implicit coupling to calculate the particle-fluid force.
   *
   * @param[in] void_fraction_solution The void fraction value calculated with
   * one of the methods of the VoidFractionBase class.
   *
   * @param[in] particle_handler The particle handler object that stores and
   * manages the particles in the simulations.
   *
   * @param[in] current_filtered_VOF_solution The present value of the VOF
   * solution.
   *
   * @param[in] drag_coupling Indicator for the type of coupling that is
   * used. This parameter essentially decides which velocity_pressure solution
   * is interpolated at the location of the particles.
   */
  template <typename VectorType>
  void
  reinit_particle_fluid_interactions(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
    const typename DoFHandler<dim>::active_cell_iterator &phase_cell,
    const VectorType                      &present_velocity_pressure_solution,
    const VectorType                      &previous_velocity_pressure_solution,
    const VectorType                      &void_fraction_solution,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::DragCoupling        &drag_coupling,
    const VectorType                      &current_filtered_VOF_solution)
  {
    pic = particle_handler.particles_in_cell(velocity_cell);

    extract_particle_properties();
    reinit_particle_fluid_forces();
    calculate_cell_void_fraction(total_particle_volume);

    if (number_of_particles == 0)
      return;

    Quadrature<dim> q_particles_location =
      gather_particles_reference_location();
    calculate_fluid_fields_at_particle_location(
      q_particles_location,
      velocity_cell,
      present_velocity_pressure_solution,
      previous_velocity_pressure_solution,
      drag_coupling);

    if (this->interpolated_void_fraction)
      {
        calculate_void_fraction_at_particle_location(q_particles_location,
                                                     void_fraction_cell,
                                                     void_fraction_solution);
      }
    calculate_vof_at_particle_location(q_particles_location,
                                       phase_cell,
                                       current_filtered_VOF_solution);
    calculate_fluid_properties_at_particle_location();
    calculate_force_parameters_at_particle_location();
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

  /**
   * @brief Reinitialize the content of the scratch for the heat transfer
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
    if (time_stepping_is_bdf(method))
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

  /**
   * @brief Reinitialize the content of the scratch for CH
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the CH FE and not the
   * Navier-Stokes FE
   *
   * @param current_solution The present value of the solution for [phi]
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

    // Gather filtered VOF solution (values, gradients)
    this->fe_values_cahn_hilliard->operator[](phase_order)
      .get_function_values(current_filtered_solution,
                           this->filtered_phase_order_cahn_hilliard_values);
  }

  /**
   * @brief enable_mortar Enables the calculation of the rotor rotation angle
   */
  void
  enable_mortar();


  /**
   * @brief Renitialize the content of the scratch data for mortar
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the FE which is used to fill the
   * FeValues.
   *
   * @param[in] mortar_parameters Parameters for the mortar method
   *
   * @param[in] radius Radius of the rotor domain
   */
  void
  reinit_mortar(const typename DoFHandler<dim>::active_cell_iterator &cell,
                const Parameters::Mortar<dim> &mortar_parameters,
                const double                  &radius);

  /**
   * @brief Calculates the physical properties. This function calculates the
   * physical properties that may be required by the fluid dynamics problem.
   * Namely the kinematic viscosity and, when required, the density.
   */
  void
  calculate_physical_properties();

  /**
   * @brief Calculates the particle forces and velocities that were projected on
   * the fluid dofs at the quadrature points of the velocity and pressure FE.
   * The values are stored in the corresponding variables in scratch data.
   *
   * @param[in] particle_drag_cell Iterator pointing to the current active cell
   * using the particle drag force DoFHandler
   *
   * @param[in] particle_two_way_coupling_force_cell Iterator pointing to the
   * current active cell using the particle two-way coupling force DoFHandler
   *
   * @param[in] particle_velocity_cell Iterator pointing to the current active
   * cell using the particle velocity DoFHandler
   *
   * @param[in] particle_momentum_transfer_coefficient_cell Iterator pointing to
   * the current active cell using the particle momentum transfer coefficient
   * DoFHandler
   *
   * @param[in] particle_fluid_drag Object containing the projection of the drag
   * calculated for the particles onto the fluid dofs
   *
   * @param[in] particle_fluid_force_two_way_coupling Object containing the
   * projection of the two-way coupling force calculated for the particles onto
   * the fluid dofs
   *
   * @param[in] particle_velocity Object containing the projection of the
   * particle velocities onto the fluid dofs
   *
   * @param[in] particle_momentum_transfer_coefficient Object containing the
   * projection of the momentum transfer coefficient calculated for the
   * particles onto the fluid dofs
   *
   * @param[in] drag_coupling Enumeration specifying the numerical coupling
   * strategy used for the computation of the drag force between the fluid and
   * the particles.
   */
  template <typename VectorType>
  void
  calculate_particle_fields_values(
    const typename DoFHandler<dim>::active_cell_iterator &particle_drag_cell,
    const typename DoFHandler<dim>::active_cell_iterator
      &particle_two_way_coupling_force_cell,
    const typename DoFHandler<dim>::active_cell_iterator
      &particle_velocity_cell,
    const typename DoFHandler<dim>::active_cell_iterator
                                   &particle_momentum_transfer_coefficient_cell,
    const VectorType               &particle_fluid_drag,
    const VectorType               &particle_fluid_force_two_way_coupling,
    const VectorType               &particle_velocity,
    const VectorType               &particle_momentum_transfer_coefficient,
    const Parameters::DragCoupling &drag_coupling)
  {
    constexpr FEValuesExtractors::Vector vector_index(0);

    this->fe_values_particle_two_way_coupling_force->reinit(
      particle_two_way_coupling_force_cell);
    (*this->fe_values_particle_two_way_coupling_force)[vector_index]
      .get_function_values(particle_fluid_force_two_way_coupling,
                           this->particle_two_way_coupling_force_values);

    // particle_drag_values will remain zero in the implicit and semi-implicit
    // coupling, since the momentum transfer coefficient is used instead, while
    // particle_momentum_transfer_coefficient_values and
    // particle_velocity_values will remain zero in the fully explicit coupling
    if (drag_coupling == Parameters::DragCoupling::fully_explicit)
      {
        this->fe_values_particle_drag->reinit(particle_drag_cell);
        (*this->fe_values_particle_drag)[vector_index].get_function_values(
          particle_fluid_drag, this->particle_drag_values);
      }
    else
      {
        this->fe_values_particle_velocity->reinit(particle_velocity_cell);
        (*this->fe_values_particle_velocity)[vector_index].get_function_values(
          particle_velocity, this->particle_velocity_values);

        this->fe_values_particle_momentum_transfer_coefficient->reinit(
          particle_momentum_transfer_coefficient_cell);
        (*this->fe_values_particle_momentum_transfer_coefficient)
          .get_function_values(
            particle_momentum_transfer_coefficient,
            this->particle_momentum_transfer_coefficient_values);
      }
  }

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
  double              density_scale;
  // The density at the particles locations is used to calculate the forces
  // between the fluid and the particles.
  std::vector<double> density_at_particle_location;
  std::vector<double> compressibility_multiplier;
  std::vector<double> dynamic_viscosity_0;
  std::vector<double> dynamic_viscosity_1;
  std::vector<double> kinematic_viscosity_0;
  std::vector<double> kinematic_viscosity_1;
  /// Scale of the kinematic viscosity for fluid 0.
  double kinematic_viscosity_scale_0;
  /// Scale of the kinematic viscosity for fluid 1.
  double kinematic_viscosity_scale_1;
  /// Kinematic viscosity at the particles locations. This is used to calculate
  /// the Reynolds number and the forces between the fluid and the particles.
  std::vector<double> kinematic_viscosity_at_particle_location;
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
  Table<1, double>        JxW;
  std::vector<Point<dim>> quadrature_points;

  // Components index
  std::vector<unsigned int> components;

  // Velocity and pressure values
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<double>                      velocity_divergences;
  std::vector<Tensor<2, dim>>              velocity_gradients;
  std::vector<Tensor<1, dim>>              velocity_laplacians;
  std::vector<Tensor<3, dim>>              velocity_hessians;
  std::vector<Tensor<1, dim>>              velocity_for_stabilization;
  std::vector<double>                      shear_rate;
  std::vector<double>                      pressure_values;
  std::vector<Tensor<1, dim>>              pressure_gradients;
  std::vector<std::vector<double>>         previous_pressure_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<Tensor<1, dim>>              sdirk_stage_sum;

  // Shape functions
  Table<2, double>         div_phi_u;
  Table<2, Tensor<1, dim>> phi_u;
  Table<2, Tensor<3, dim>> hess_phi_u;
  Table<2, Tensor<1, dim>> laplacian_phi_u;
  Table<2, Tensor<2, dim>> grad_phi_u;
  Table<2, double>         phi_p;
  Table<2, Tensor<1, dim>> grad_phi_p;

  /**
   * Scratch component for the VOF auxiliary physics
   */
  bool                             gather_vof;
  unsigned int                     n_dofs_vof;
  std::vector<double>              phase_values;
  std::vector<double>              filtered_phase_values;
  std::vector<double>              filtered_phase_values_at_particle_location;
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
  bool gather_particles_information;
  bool interpolated_void_fraction;

  std::shared_ptr<FEValues<dim>> fe_values_particle_drag;
  std::shared_ptr<FEValues<dim>> fe_values_particle_two_way_coupling_force;
  std::shared_ptr<FEValues<dim>> fe_values_particle_velocity;
  std::shared_ptr<FEValues<dim>>
    fe_values_particle_momentum_transfer_coefficient;

  bool                        gather_particle_field_project;
  std::vector<Tensor<1, dim>> particle_velocity;
  std::vector<Tensor<1, dim>> particle_velocity_values;
  std::vector<Tensor<1, dim>> particle_drag_values;
  std::vector<Tensor<1, dim>> particle_two_way_coupling_force_values;
  std::vector<double>         particle_momentum_transfer_coefficient_values;
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
  double         total_particle_volume;
  double         cell_volume;
  double         beta_drag;
  Tensor<1, dim> explicit_particle_volumetric_acceleration_on_fluid;


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
   * Scratch component for the mortar method
   */
  bool                        gather_mortar;
  std::vector<Tensor<1, dim>> rotor_linear_velocity_values;

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
  Table<2, double>                     face_JxW;
  std::vector<std::vector<Point<dim>>> face_quadrature_points;
  Table<2, Tensor<1, dim>>             face_normal;

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
  Table<3, double>         face_div_phi_u;
  Table<3, Tensor<1, dim>> face_phi_u;
  Table<3, Tensor<3, dim>> face_hess_phi_u;
  Table<3, Tensor<1, dim>> face_laplacian_phi_u;
  Table<3, Tensor<2, dim>> face_grad_phi_u;
  Table<3, double>         face_phi_p;
  Table<3, Tensor<1, dim>> face_grad_phi_p;

  /// SDIRK Butcher table with coefficients for the SDIRK time-stepping method
  SDIRKTable sdirk_table;
};

#endif
