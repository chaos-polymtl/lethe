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
 * Scratch data for the CahnHilliard auxiliary physics
 */

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


#ifndef lethe_cahn_hilliard_scratch_data_h
#  define lethe_cahn_hilliard_scratch_data_h

using namespace dealii;


/**
 * @brief Class that stores the information required by the assembly procedure
 * for the Cahn-Hilliard equations. Consequently, this class
 * calculates the phase field parameter Phi (values, gradients, laplacians),
 * the chemical potential eta (values, gradients, laplacians) and the shape
 * function (values, gradients, laplacians) at all the Gauss points for all
 * degrees of freedom and stores it into arrays. Additionally, the user can
 * request that this class gathers additional fields for physics which are
 * coupled to the Cahn-Hilliard equations, such as the velocity which is
 * required. This class serves as a separation between the evaluation at the
 * Gauss point of the variables of interest and their use in the assembly, which
 * is carried out by the assembler functions. For more information on this
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
class CahnHilliardScratchData
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
   * @param fe_cahn_hilliard The FESystem used to solve the Cahn-Hilliard equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Cahn-Hilliard equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  CahnHilliardScratchData(const PhysicalPropertiesManager &properties_manager,
                          const FESystem<dim>             &fe_cahn_hilliard,
                          const Quadrature<dim>           &quadrature,
                          const Mapping<dim>              &mapping,
                          const FiniteElement<dim>        &fe_fd,
                          const Quadrature<dim - 1>       &face_quadrature)
    : properties_manager(properties_manager)
    , fe_values_cahn_hilliard(mapping,
                              fe_cahn_hilliard,
                              quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values | update_gradients |
                                update_hessians)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values)
    , fe_face_values_cahn_hilliard(mapping,
                                   fe_cahn_hilliard,
                                   face_quadrature,
                                   update_values | update_quadrature_points |
                                     update_JxW_values | update_gradients |
                                     update_normal_vectors)
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
  CahnHilliardScratchData(const CahnHilliardScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_cahn_hilliard(sd.fe_values_cahn_hilliard.get_mapping(),
                              sd.fe_values_cahn_hilliard.get_fe(),
                              sd.fe_values_cahn_hilliard.get_quadrature(),
                              update_values | update_quadrature_points |
                                update_JxW_values | update_gradients |
                                update_hessians)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values)
    , fe_face_values_cahn_hilliard(
        sd.fe_face_values_cahn_hilliard.get_mapping(),
        sd.fe_face_values_cahn_hilliard.get_fe(),
        sd.fe_face_values_cahn_hilliard.get_quadrature(),
        update_values | update_quadrature_points | update_JxW_values |
          update_gradients | update_normal_vectors)
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
   * @param source_function The function describing the source term in Cahn-Hilliard
   * equations
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         Function<dim>                 *source_function,
         Parameters::CahnHilliard       cahn_hilliard_parameters)
  {
    this->phase_order.component        = 0;
    this->chemical_potential.component = 1;

    this->fe_values_cahn_hilliard.reinit(cell);

    quadrature_points = this->fe_values_cahn_hilliard.get_quadrature_points();
    auto &fe_cahn_hilliard = this->fe_values_cahn_hilliard.get_fe();

    source_function->value_list(quadrature_points, source_phase_order, 0);
    source_function->value_list(quadrature_points,
                                source_chemical_potential,
                                1);


    if (dim == 2)
      this->cell_size =
        std::sqrt(4. * cell->measure() / M_PI) / fe_cahn_hilliard.degree;
    else if (dim == 3)
      this->cell_size =
        pow(6 * cell->measure() / M_PI, 1. / 3.) / fe_cahn_hilliard.degree;

    // Gather Phi and eta (values, gradient and laplacian)
    this->fe_values_cahn_hilliard[phase_order].get_function_values(
      current_solution, this->phase_order_values);
    this->fe_values_cahn_hilliard[phase_order].get_function_gradients(
      current_solution, this->phase_order_gradients);
    this->fe_values_cahn_hilliard[phase_order].get_function_laplacians(
      current_solution, this->phase_order_laplacians);

    this->fe_values_cahn_hilliard[chemical_potential].get_function_values(
      current_solution, this->chemical_potential_values);
    this->fe_values_cahn_hilliard[chemical_potential].get_function_gradients(
      current_solution, this->chemical_potential_gradients);
    this->fe_values_cahn_hilliard[chemical_potential].get_function_laplacians(
      current_solution, this->chemical_potential_laplacians);


    // Gather previous phase order values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_cahn_hilliard[phase_order].get_function_values(
          previous_solutions[p], previous_phase_order_values[p]);
      }

    // Gather previous chemical potential values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_cahn_hilliard[chemical_potential].get_function_values(
          previous_solutions[p], previous_chemical_potential_values[p]);
      }

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_cahn_hilliard.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Shape functions for the phase order
            this->phi_phase[q][k] =
              this->fe_values_cahn_hilliard[phase_order].value(k, q);
            this->grad_phi_phase[q][k] =
              this->fe_values_cahn_hilliard[phase_order].gradient(k, q);
            this->hess_phi_phase[q][k] =
              this->fe_values_cahn_hilliard[phase_order].hessian(k, q);
            this->laplacian_phi_phase[q][k] = trace(this->hess_phi_phase[q][k]);

            // Shape functions for the chemical potential
            this->phi_potential[q][k] =
              this->fe_values_cahn_hilliard[chemical_potential].value(k, q);
            this->grad_phi_potential[q][k] =
              this->fe_values_cahn_hilliard[chemical_potential].gradient(k, q);
            this->hess_phi_potential[q][k] =
              this->fe_values_cahn_hilliard[chemical_potential].hessian(k, q);
            this->laplacian_phi_potential[q][k] =
              trace(this->hess_phi_potential[q][k]);
          }
      }

    this->is_boundary_cell =
      cell->at_boundary(); // The attribute needs to be updated because the
                           // assembler for the angle of contact boundary
                           // condition needs to know if the cell is at the
                           // boundary
    if (this->is_boundary_cell)
      {
        n_faces          = cell->n_faces();
        is_boundary_face = std::vector<bool>(n_faces, false);
        n_faces_q_points = fe_face_values_cahn_hilliard.get_quadrature().size();
        boundary_face_id = std::vector<unsigned int>(n_faces);

        face_JxW = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        this->grad_phi_face_phase =
          std::vector<std::vector<std::vector<Tensor<1, dim>>>>(
            n_faces,
            std::vector<std::vector<Tensor<1, dim>>>(
              n_faces_q_points, std::vector<Tensor<1, dim>>(n_dofs)));

        this->face_phase_grad_values = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        this->face_normal = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

        for (const auto face : cell->face_indices())
          {
            this->is_boundary_face[face] = cell->face(face)->at_boundary();
            if (this->is_boundary_face[face])
              {
                fe_face_values_cahn_hilliard.reinit(cell, face);
                boundary_face_id[face] = cell->face(face)->boundary_id();

                this->fe_face_values_cahn_hilliard[phase_order]
                  .get_function_gradients(current_solution,
                                          this->face_phase_grad_values[face]);

                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  {
                    face_JxW[face][q] = fe_face_values_cahn_hilliard.JxW(q);
                    this->face_normal[face][q] =
                      this->fe_face_values_cahn_hilliard.normal_vector(q);
                    for (const unsigned int k :
                         fe_face_values_cahn_hilliard.dof_indices())
                      {
                        this->grad_phi_face_phase[face][q][k] =
                          this->fe_face_values_cahn_hilliard[phase_order]
                            .gradient(k, q);
                      }
                  }
              }
          }
      }

    // CH epsilon parameter
    this->epsilon = (cahn_hilliard_parameters.epsilon_set_method ==
                     Parameters::EpsilonSetStrategy::manual) ?
                      cahn_hilliard_parameters.epsilon :
                      2 * this->cell_size;
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


  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  dealii::types::material_id           material_id;
  double                               epsilon;
  std::vector<double>                  density;
  std::vector<double>                  kinematic_viscosity;


  FEValuesExtractors::Scalar phase_order;
  FEValuesExtractors::Scalar chemical_potential;

  // FEValues for the Cahn-Hilliard problem
  FEValues<dim> fe_values_cahn_hilliard;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Phase order and chemical potential values
  std::vector<double>              phase_order_values;
  std::vector<Tensor<1, dim>>      phase_order_gradients;
  std::vector<double>              phase_order_laplacians;
  std::vector<std::vector<double>> previous_phase_order_values;
  std::vector<double>              chemical_potential_values;
  std::vector<Tensor<1, dim>>      chemical_potential_gradients;
  std::vector<double>              chemical_potential_laplacians;
  std::vector<std::vector<double>> previous_chemical_potential_values;

  // Source term
  std::vector<double> source_phase_order;
  std::vector<double> source_chemical_potential;

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
  // This FEValues must be instantiated for the velocity
  FEValues<dim>                            fe_values_fd;
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<Tensor<2, dim>>              velocity_gradient_values;

  // Scratch for the face boundary condition
  FEFaceValues<dim>                fe_face_values_cahn_hilliard;
  std::vector<std::vector<double>> face_JxW;

  unsigned int n_faces;
  unsigned int n_faces_q_points;

  // Is boundary cell indicator
  bool                      is_boundary_cell;
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;

  // First vector is face number, second quadrature point, third DOF
  std::vector<std::vector<std::vector<Tensor<1, dim>>>> grad_phi_face_phase;
  // First vector is face number, second quadrature point
  std::vector<std::vector<Tensor<1, dim>>> face_phase_grad_values;
  // The normal vector is necessary for the free angle boundary condition
  std::vector<std::vector<Tensor<1, dim>>> face_normal;
};

#endif
