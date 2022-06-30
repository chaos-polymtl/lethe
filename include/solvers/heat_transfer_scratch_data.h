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
 * Implementation of heat transfer as an auxiliary physics.
 * This heat equation is weakly coupled to the velocity field.
 * Equation solved:
 * rho * Cp * (dT/dt + u.gradT) = k div(gradT) + nu/rho * (gradu : gradu)
 *
 * Polytechnique Montreal, 2020-
 */

#ifndef lethe_heat_transfer_scratch_data_h
#define lethe_heat_transfer_scratch_data_h

#include <core/density_model.h>
#include <core/multiphysics.h>
#include <core/physical_property_model.h>
#include <core/specific_heat_model.h>
#include <core/thermal_conductivity_model.h>

#include <solvers/multiphysics_interface.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/numerics/vector_tools.h>

using namespace dealii;


/**
 * @brief HeatTransferScratchData class
 * stores the information required by the assembly procedure
 * for a heat transfer advection-diffusion equation. Consequently, this class
 * calculates the heat transfer (values, gradients, laplacians) and the shape
 * function (values, gradients, laplacians) at all the gauss points for all
 * degrees of freedom and stores it into arrays. Additionnaly, the user can
 * request that this class gathers additional fields for physics which are
 * coupled to the heat transfer equation, such as the velocity which is
 * required. This class serves as a seperation between the evaluation at the
 * gauss point of the variables of interest and their use in the assembly, which
 * is carried out by the assembler methods. For more information on this
 * design, the reader can consult deal.II step-9
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
class HeatTransferScratchData
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evalution, since this needs to be done at the cell level.
   *
   * @param properties_manager The physical properties Manager (see physical_properties_manager.h)
   *
   * @param fe_ht The FESystem used to solve the Heat Transfer equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   * @param face_quadrature_formula The face quadrature formula
   *
   */
  HeatTransferScratchData(const PhysicalPropertiesManager properties_manager,
                          const FiniteElement<dim> &      fe_ht,
                          const Quadrature<dim> &         quadrature,
                          const Mapping<dim> &            mapping,
                          const FiniteElement<dim> &      fe_fd,
                          const Quadrature<dim - 1> &     face_quadrature)
    : properties_manager(properties_manager)
    , fe_values_T(mapping,
                  fe_ht,
                  quadrature,
                  update_values | update_quadrature_points | update_JxW_values |
                    update_gradients | update_hessians)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
    , fe_face_values_ht(mapping,
                        fe_ht,
                        face_quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values)
  {
    gather_vof = false;

    allocate();
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism, it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd The scratch data
   */
  HeatTransferScratchData(const HeatTransferScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_T(sd.fe_values_T.get_mapping(),
                  sd.fe_values_T.get_fe(),
                  sd.fe_values_T.get_quadrature(),
                  update_values | update_quadrature_points | update_JxW_values |
                    update_gradients | update_hessians)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values | update_gradients)
    , fe_face_values_ht(sd.fe_face_values_ht.get_mapping(),
                        sd.fe_face_values_ht.get_fe(),
                        sd.fe_face_values_ht.get_quadrature(),
                        update_values | update_quadrature_points |
                          update_JxW_values)
  {
    gather_vof = sd.gather_vof;
    allocate();
    if (sd.gather_vof)
      enable_vof(sd.fe_values_vof->get_fe(),
                 sd.fe_values_vof->get_quadrature(),
                 sd.fe_values_vof->get_mapping());
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
   * Using the FeValues and the content of the solutions, previous solutions and
   * solutions stages, fills all of the class member of the scratch
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the FeValues
   *
   * @param current_solution The present value of the solution for the Heat Transfer
   *
   * @param previous_solutions The solutions at the previous time steps
   *
   * @param solution_stages The solution at the intermediary stages (for SDIRK methods)
   *
   * @param source_function The function describing the Heat Transfer source term
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType &                                    current_solution,
         const std::vector<TrilinosWrappers::MPI::Vector> &previous_solutions,
         const std::vector<VectorType> &                   solution_stages,
         Function<dim> *                                   source_function)
  {
    this->fe_values_T.reinit(cell);

    quadrature_points = this->fe_values_T.get_quadrature_points();
    auto &fe_T        = this->fe_values_T.get_fe();

    source_function->value_list(quadrature_points, source);

    if (dim == 2)
      this->cell_size = std::sqrt(4. * cell->measure() / M_PI) / fe_T.degree;
    else if (dim == 3)
      this->cell_size = pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_T.degree;

    // Gather temperature (values, gradient and laplacian)
    this->fe_values_T.get_function_values(current_solution,
                                          this->present_temperature_values);
    this->fe_values_T.get_function_gradients(current_solution,
                                             this->temperature_gradients);
    this->fe_values_T.get_function_laplacians(
      current_solution, this->present_temperature_laplacians);

    // Gather previous temperature values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_T.get_function_values(previous_solutions[p],
                                              previous_temperature_values[p]);

        this->fe_values_T.get_function_gradients(
          previous_solutions[p], previous_temperature_gradients[p]);
      }

    // Gather temperature stages
    for (unsigned int s = 0; s < solution_stages.size(); ++s)
      {
        this->fe_values_T.get_function_values(solution_stages[s],
                                              stages_temperature_values[s]);
      }

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_T.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Shape function
            this->phi_T[q][k]           = this->fe_values_T.shape_value(k, q);
            this->grad_phi_T[q][k]      = this->fe_values_T.shape_grad(k, q);
            this->hess_phi_T[q][k]      = this->fe_values_T.shape_hessian(k, q);
            this->laplacian_phi_T[q][k] = trace(this->hess_phi_T[q][k]);
          }
      }


    // Arrays related to faces must be re-initialized for each cell, since they
    // might depend on reference cell
    // Only carry out this initialization if the cell is a boundary cell,
    // otherwise these are wasted calculations
    this->is_boundary_cell = cell->at_boundary();
    if (cell->at_boundary())
      {
        n_faces          = cell->n_faces();
        is_boundary_face = std::vector<bool>(n_faces, false);
        n_faces_q_points = fe_face_values_ht.get_quadrature().size();
        boundary_face_id = std::vector<unsigned int>(n_faces);

        face_JxW = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));


        this->phi_face_T = std::vector<std::vector<std::vector<double>>>(
          n_faces,
          std::vector<std::vector<double>>(n_faces_q_points,
                                           std::vector<double>(n_dofs)));

        this->temperature_face_value = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        for (const auto face : cell->face_indices())
          {
            this->is_boundary_face[face] = cell->face(face)->at_boundary();
            if (this->is_boundary_face[face])
              {
                fe_face_values_ht.reinit(cell, face);
                boundary_face_id[face] = cell->face(face)->boundary_id();
                this->fe_face_values_ht.get_function_values(
                  current_solution, this->temperature_face_value[face]);

                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  {
                    face_JxW[face][q] = fe_face_values_ht.JxW(q);
                    for (const unsigned int k : fe_face_values_ht.dof_indices())
                      {
                        this->phi_face_T[face][q][k] =
                          this->fe_face_values_ht.shape_value(k, q);
                      }
                  }
              }
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
                  const VectorType &current_solution)
  {
    this->fe_values_fd.reinit(cell);

    this->fe_values_fd[velocities].get_function_values(current_solution,
                                                       velocity_values);
  }

  /** @brief Reinitialize the velocity gradient, calculated by the Fluid Dynamics
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param current_solution The present value of the solution for [u,p]
   *
   */

  template <typename VectorType>
  void
  reinit_velocity_gradient(const VectorType &current_solution)
  {
    this->fe_values_fd[velocities].get_function_gradients(
      current_solution, velocity_gradient_values);
  }


  /**
   * @brief enable_vof Enables the collection of the VOF data by the scratch.
   *
   * @param fe FiniteElement associated with the VOF.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly.
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly.
   */

  void
  enable_vof(const FiniteElement<dim> &fe,
             const Quadrature<dim> &   quadrature,
             const Mapping<dim> &      mapping);

  /** @brief Reinitialize the content of the scratch for VOF.
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the VOF FE and not the
   * Fluid Dynamics FE
   *
   * @param current_solution The present solution for the phase value
   *
   * @param solution_stages The solution at the intermediary stages
   * (for SDIRK methods) for the phase value
   *
   * NB: the previous_solutions are not used in heat_transfer_assemblers,
   * contrary to navier_stokes_assemblers
   */

  template <typename VectorType>
  void
  reinit_vof(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType &current_solution,
             const std::vector<VectorType> & /*solution_stages*/)
  {
    this->fe_values_vof->reinit(cell);
    // Gather phase fraction (values, gradient)
    this->fe_values_vof->get_function_values(current_solution,
                                             this->phase_values);
  }


  /** @brief Calculates the physical properties. This method calculates the physical properties
   * that may be required by the heat transfer problem. Namely the density,
   * specific heat, thermal conductivity and viscosity (for viscous
   * dissipation).
   *
   */
  void
  calculate_physical_properties();


  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  std::vector<double>                  specific_heat;
  std::vector<double>                  thermal_conductivity;
  std::vector<double>                  density;
  std::vector<double>                  viscosity;
  // Gradient of the specific heat with respect to the temperature
  // This is calculated by deriving the specific heat by the temperature
  // (dCp/dT)
  std::vector<double> grad_specific_heat_temperature;

  // Auxiliary property vector for VOF simulations
  std::vector<double> specific_heat_0;
  std::vector<double> thermal_conductivity_0;
  std::vector<double> density_0;
  std::vector<double> viscosity_0;
  std::vector<double> grad_specific_heat_temperature_0;

  std::vector<double> specific_heat_1;
  std::vector<double> thermal_conductivity_1;
  std::vector<double> density_1;
  std::vector<double> viscosity_1;
  std::vector<double> grad_specific_heat_temperature_1;


  // FEValues for the HT problem
  FEValues<dim> fe_values_T;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;



  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Temperature values
  std::vector<double>                      present_temperature_values;
  std::vector<Tensor<1, dim>>              temperature_gradients;
  std::vector<double>                      present_temperature_laplacians;
  std::vector<double>                      present_face_temperature_values;
  std::vector<std::vector<double>>         previous_temperature_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_temperature_gradients;
  std::vector<std::vector<double>>         stages_temperature_values;


  // Shape functions and gradients
  std::vector<std::vector<double>>         phi_T;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi_T;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi_T;
  std::vector<std::vector<double>>         laplacian_phi_T;

  // Source term
  std::vector<double> source;

  /**
   * Scratch component for the VOF auxiliary physics
   */
  bool                gather_vof;
  unsigned int        n_dofs_vof;
  std::vector<double> phase_values;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_vof;

  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<Tensor<2, dim>> velocity_gradient_values;
  std::vector<double>         shear_rate_values;

  // Scratch for the face boundary condition
  FEFaceValues<dim>                fe_face_values_ht;
  std::vector<std::vector<double>> face_JxW;

  unsigned int n_faces;
  unsigned int n_faces_q_points;

  // Is boundary cell indicator
  bool                      is_boundary_cell;
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;


  // First vector is face number, second quadrature point, third DOF
  std::vector<std::vector<std::vector<double>>> phi_face_T;
  // First vector is face number, second quadrature point
  std::vector<std::vector<double>> temperature_face_value;
};

#endif
