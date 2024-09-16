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
 */


#ifndef lethe_tracer_scratch_data_h
#define lethe_tracer_scratch_data_h

#include <core/ale.h>
#include <core/multiphysics.h>

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
 * @brief Class that stores the information required by the assembly procedure
 * for a Tracer advection-diffusion equation. Consequently, this class
 *calculates the tracer (values, gradients, laplacians) and the shape function
 * (values, gradients, laplacians) at all the gauss points for all degrees
 * of freedom and stores it into arrays. Additionnaly, the use can request
 * that this class gathers additional fields for physics which are coupled
 * to the Tracer equation, such as the velocity which is required. This class
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
class TracerScratchData
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
  TracerScratchData(const PhysicalPropertiesManager &properties_manager,
                    const FiniteElement<dim>        &fe_tracer,
                    const Quadrature<dim>           &quadrature,
                    const Mapping<dim>              &mapping,
                    const FiniteElement<dim>        &fe_fd,
                    const Quadrature<dim - 1>       &face_quadrature)
    : properties_manager(properties_manager)
    , fe_values_tracer(mapping,
                       fe_tracer,
                       quadrature,
                       update_values | update_quadrature_points |
                         update_JxW_values | update_gradients | update_hessians)
    , fe_face_values_tracer(mapping,
                            fe_tracer,
                            face_quadrature,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_normal_vectors)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values)
    , fe_face_values_fd(mapping,
                        fe_fd,
                        face_quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values | update_normal_vectors)
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
  TracerScratchData(const TracerScratchData<dim> &sd)
    : properties_manager(sd.properties_manager)
    , fe_values_tracer(sd.fe_values_tracer.get_mapping(),
                       sd.fe_values_tracer.get_fe(),
                       sd.fe_values_tracer.get_quadrature(),
                       update_values | update_quadrature_points |
                         update_JxW_values | update_gradients | update_hessians)
    , fe_face_values_tracer(sd.fe_face_values_tracer.get_mapping(),
                            sd.fe_face_values_tracer.get_fe(),
                            sd.fe_face_values_tracer.get_quadrature(),
                            update_values | update_quadrature_points |
                              update_JxW_values | update_normal_vectors)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values)
    , fe_face_values_fd(sd.fe_face_values_fd.get_mapping(),
                        sd.fe_face_values_fd.get_fe(),
                        sd.fe_face_values_fd.get_quadrature(),
                        update_values | update_quadrature_points |
                          update_JxW_values | update_normal_vectors)
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

  /** @brief Reinitializes the content of the scratch.
   *
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the
   * FeValues.
   *
   * @param[in] current_solution The present value of the solution.
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   * @param[in] source_function The function describing the tracer source term.
   *
   * @param[in] levelset_function The function describing the particles (if
   * there are any).
   */
  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         const Function<dim>           *source_function,
         const Function<dim>           *levelset_function)
  {
    this->fe_values_tracer.reinit(cell);

    quadrature_points = this->fe_values_tracer.get_quadrature_points();
    auto &fe_tracer   = this->fe_values_tracer.get_fe();

    source_function->value_list(quadrature_points, source);

    if (properties_manager.field_is_required(field::levelset))
      {
        Assert(
          levelset_function != nullptr,
          ExcMessage(
            "Levelset function is required for tracer assembly, but the level set function is a nullptr"));

        levelset_function->value_list(quadrature_points, sdf_values);
      }

    // Compute cell diameter
    double cell_measure =
      compute_cell_measure_with_JxW(this->fe_values_tracer.get_JxW_values());
    this->cell_size =
      compute_cell_diameter<dim>(cell_measure, fe_tracer.degree);

    // Gather tracer (values, gradient and laplacian)
    this->fe_values_tracer.get_function_values(current_solution,
                                               this->tracer_values);
    this->fe_values_tracer.get_function_gradients(current_solution,
                                                  this->tracer_gradients);
    this->fe_values_tracer.get_function_laplacians(current_solution,
                                                   this->tracer_laplacians);

    // Gather previous tracer values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values_tracer.get_function_values(previous_solutions[p],
                                                   previous_tracer_values[p]);
      }

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        this->JxW[q] = this->fe_values_tracer.JxW(q);

        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            // Shape function
            this->phi[q][k]      = this->fe_values_tracer.shape_value(k, q);
            this->grad_phi[q][k] = this->fe_values_tracer.shape_grad(k, q);
            this->hess_phi[q][k] = this->fe_values_tracer.shape_hessian(k, q);
            this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
          }
      }


    // Arrays related to faces must be re-initialized for each cell, since they
    // might depend on reference cell
    // Only carry out this initialization if the cell is a boundary cell,
    // otherwise these are wasted calculations
    this->is_boundary_cell = cell->at_boundary();
    if (this->is_boundary_cell)
      {
        n_faces          = cell->n_faces();
        is_boundary_face = std::vector<bool>(n_faces, false);
        n_faces_q_points = fe_face_values_tracer.get_quadrature().size();
        boundary_face_id = std::vector<unsigned int>(n_faces);

        face_JxW = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        this->phi_face_tracer = std::vector<std::vector<std::vector<double>>>(
          n_faces,
          std::vector<std::vector<double>>(n_faces_q_points,
                                           std::vector<double>(face_n_dofs)));

        this->tracer_face_value = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        for (const auto face : cell->face_indices())
          {
            this->is_boundary_face[face] = cell->face(face)->at_boundary();
            boundary_face_id[face]       = cell->face(face)->boundary_id();
            if (this->is_boundary_face[face])
              {
                fe_face_values_tracer.reinit(cell, face);
                this->fe_face_values_tracer.get_function_values(
                  current_solution, this->tracer_face_value[face]);
                face_JxW[face] = fe_face_values_tracer.get_JxW_values();
                for (unsigned int q = 0; q < n_faces_q_points; ++q)
                  {
                    for (const unsigned int k :
                         fe_face_values_tracer.dof_indices())
                      {
                        this->phi_face_tracer[face][q][k] =
                          this->fe_face_values_tracer.shape_value(k, q);
                      }
                  }
              }
          }
      }
  }

  /** @brief Reinitialize the velocity, calculated by the fluid dynamics while also taking into account ALE
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell for which the velocity is reinitialized
   * This cell must be compatible with the Fluid Dynamics FE
   *
   * @param current_solution The present value of the solution for [u,p]
   *
   * @param ale The ALE parameters which include the ALE function
   *
   * @param drift_velocity Function to calculate the drift velocity
   *
   */

  template <typename VectorType>
  void
  reinit_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const Parameters::ALE<dim>                           &ale,
    std::shared_ptr<Functions::ParsedFunction<dim>>       drift_velocity)
  {
    this->fe_values_fd.reinit(cell);

    this->fe_values_fd[velocities].get_function_values(current_solution,
                                                       velocity_values);

    // Add the drift velocity to the velocity to account for tracer drift flux
    // modeling
    Tensor<1, dim> drift_velocity_tensor;
    Vector<double> drift_velocity_vector(dim);

    for (unsigned int q = 0; q < n_q_points; ++q)
      {
        drift_velocity->vector_value(quadrature_points[q],
                                     drift_velocity_vector);
        for (unsigned int d = 0; d < dim; ++d)
          drift_velocity_tensor[d] = drift_velocity_vector[d];

        velocity_values[q] += drift_velocity_tensor;
      }

    if (ale.enabled())
      {
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

    is_boundary_cell = cell->at_boundary();
    if (is_boundary_cell)
      {
        n_faces          = cell->n_faces();
        is_boundary_face = std::vector<bool>(n_faces, false);
        n_faces_q_points = fe_face_values_fd.get_quadrature().size();
        boundary_face_id = std::vector<unsigned int>(n_faces);

        this->face_JxW = std::vector<std::vector<double>>(
          n_faces, std::vector<double>(n_faces_q_points));

        // Velocity and pressure values
        // First vector is face number, second quadrature point
        this->face_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces,
          std::vector<Tensor<1, dim>>(n_faces_q_points, Tensor<1, dim>()));

        this->face_normal = std::vector<std::vector<Tensor<1, dim>>>(
          n_faces,
          std::vector<Tensor<1, dim>>(n_faces_q_points, Tensor<1, dim>()));

        this->face_quadrature_points = std::vector<std::vector<Point<dim>>>(
          n_faces, std::vector<Point<dim>>(n_faces_q_points, Point<dim>()));

        for (const auto face : cell->face_indices())
          {
            is_boundary_face[face] = cell->face(face)->at_boundary();
            boundary_face_id[face] = cell->face(face)->boundary_id();
            if (is_boundary_face[face])
              {
                fe_face_values_fd.reinit(cell, face);
                // Gather velocity (values)
                this->fe_face_values_fd[velocities].get_function_values(
                  current_solution, this->face_velocity_values[face]);
                this->face_JxW[face] = this->fe_face_values_fd.get_JxW_values();
                this->face_normal[face] =
                  this->fe_face_values_fd.get_normal_vectors();
                this->face_quadrature_points[face] =
                  this->fe_face_values_fd.get_quadrature_points();
              }
          }
      }
  }

  /** @brief Calculates the physical properties. This function calculates the physical properties
   * that may be required by the tracer. Namely the diffusivity.
   *
   */
  void
  calculate_physical_properties();

  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  std::vector<double>                  tracer_diffusivity;
  std::vector<double>                  tracer_diffusivity_0;
  std::vector<double>                  tracer_diffusivity_1;



  // FEValues for the Tracer problem
  FEValues<dim> fe_values_tracer;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // Tracer values
  std::vector<double>              tracer_values;
  std::vector<Tensor<1, dim>>      tracer_gradients;
  std::vector<double>              tracer_laplacians;
  std::vector<std::vector<double>> previous_tracer_values;

  // Solid signed distance function
  std::vector<double> sdf_values;

  // Source term
  std::vector<double> source;

  // Scratch for the face boundary condition
  FEFaceValues<dim>                    fe_face_values_tracer;
  std::vector<std::vector<double>>     face_JxW;
  std::vector<std::vector<Point<dim>>> face_quadrature_points;

  unsigned int n_faces;
  unsigned int n_faces_q_points;
  unsigned int face_n_dofs;

  // Is boundary cell indicator
  bool                      is_boundary_cell;
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi;
  std::vector<std::vector<double>>         laplacian_phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;

  // First vector is face number, second quadrature point, third DOF
  std::vector<std::vector<std::vector<double>>> phi_face_tracer;
  // First vector is face number, second quadrature point
  std::vector<std::vector<double>> tracer_face_value;

  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>                            fe_values_fd;
  FEFaceValues<dim>                        fe_face_values_fd;
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> face_velocity_values;
  std::vector<std::vector<Tensor<1, dim>>>
    face_normal; // TODO Initialize properly
};

#endif
