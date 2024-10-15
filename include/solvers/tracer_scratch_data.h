// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_scratch_data_h
#define lethe_tracer_scratch_data_h

#include <core/ale.h>
#include <core/multiphysics.h>

#include <solvers/physical_properties_manager.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_interface_values.h>
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
                    const Quadrature<dim - 1>       &face_quadrature,
                    const Mapping<dim>              &mapping,
                    const FiniteElement<dim>        &fe_fd)
    : properties_manager(properties_manager)
    , fe_values_tracer(mapping,
                       fe_tracer,
                       quadrature,
                       update_values | update_quadrature_points |
                         update_JxW_values | update_gradients | update_hessians)
    , fe_interface_values_tracer(mapping,
                                 fe_tracer,
                                 face_quadrature,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values | update_normal_vectors)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values)
    , fe_face_values_fd(mapping, fe_fd, face_quadrature, update_values)
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
    , fe_interface_values_tracer(sd.fe_interface_values_tracer.get_mapping(),
                                 sd.fe_interface_values_tracer.get_fe(),
                                 sd.fe_interface_values_tracer.get_quadrature(),
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values | update_normal_vectors)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values)
    , fe_face_values_fd(sd.fe_face_values_fd.get_mapping(),
                        sd.fe_face_values_fd.get_fe(),
                        sd.fe_face_values_fd.get_quadrature(),
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
    this->fe_values_tracer.get_function_gradients(
      previous_solutions[0], this->previous_tracer_gradients);

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

    boundary_index = 0;
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
    const VectorType                                     &velocity_solution,
    const Parameters::ALE<dim>                           &ale,
    std::shared_ptr<Functions::ParsedFunction<dim>>       drift_velocity)
  {
    this->fe_values_fd.reinit(cell);

    this->fe_values_fd[velocities].get_function_values(velocity_solution,
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


  /** @brief Reinitializes the content of the scratch for the internal faces. This is only used for the DG assemblers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face_no The face index associated with the cell
   *
   * @param[in] sub_face_no The subface index associated with the face
   *
   * @param[in] ncell The neighboring cell
   *
   * @param[in] n_face_no The face index associated with the neighboring cell
   *
   * @param[in] n_sub_face_no The subface index associated with the neighboring
   * cell
   *
   * @param[in] current_solution The present value of the solution.
   * (if there are any).
   *
   * @param[in] levelset_function Function used for the calculation of the
   * level set. This is used for the calculation of the physical properties.
   * there are any).
   */
  template <typename VectorType>
  void
  reinit_internal_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                   &face_no,
    const unsigned int                                   &sf,
    const typename DoFHandler<dim>::active_cell_iterator &ncell,
    const unsigned int                                   &nf,
    const unsigned int                                   &nsf,
    const VectorType                                     &current_solution,
    const Function<dim>                                  *levelset_function)
  {
    fe_interface_values_tracer.reinit(cell, face_no, sf, ncell, nf, nsf);
    face_quadrature_points = fe_interface_values_tracer.get_quadrature_points();

    n_interface_dofs = fe_interface_values_tracer.n_current_interface_dofs();

    const double extent1 = cell->measure() / cell->face(face_no)->measure();
    const double extent2 = ncell->measure() / ncell->face(nf)->measure();

    penalty_factor =
      get_penalty_factor(fe_values_tracer.get_fe().degree, extent1, extent2);


    // BB TODO : Preallocate memory here
    values_here.resize(face_quadrature_points.size());
    values_there.resize(face_quadrature_points.size());
    tracer_value_jump.resize(face_quadrature_points.size());
    tracer_average_gradient.resize(face_quadrature_points.size());


    fe_interface_values_tracer.get_fe_face_values(0).get_function_values(
      current_solution, values_here);
    fe_interface_values_tracer.get_fe_face_values(1).get_function_values(
      current_solution, values_there);

    fe_interface_values_tracer.get_jump_in_function_values(current_solution,
                                                           tracer_value_jump);

    fe_interface_values_tracer.get_average_of_function_gradients(
      current_solution, tracer_average_gradient);


    if (properties_manager.field_is_required(field::levelset))
      {
        Assert(
          levelset_function != nullptr,
          ExcMessage(
            "Levelset function is required for tracer assembly, but the level set function is a nullptr"));

        levelset_function->value_list(face_quadrature_points, sdf_face_values);
      }
  }


  /** @brief Reinitializes the content of the scratch for the internal faces. This is only used for the DG assemblers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face The face index associated with the cell
   *
   * @param[in] current_solution The present value of the solution.
   * there are any).
   *
   * @param[in] levelset_function Function used for the calculation of the level
   * set. This is used for the calculation of the physical properties. there are
   * any).
   */
  template <typename VectorType>
  void
  reinit_boundary_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                   &face,
    const unsigned int                                   &boundary_index,
    const VectorType                                     &current_solution,
    const Function<dim>                                  *levelset_function)
  {
    fe_interface_values_tracer.reinit(cell, face);
    const FEFaceValuesBase<dim> &fe_face =
      fe_interface_values_tracer.get_fe_face_values(0);
    face_quadrature_points = fe_interface_values_tracer.get_quadrature_points();

    // BB TODO: These could be pre-allocated
    values_here.resize(face_quadrature_points.size());
    gradients_here.resize(face_quadrature_points.size());

    fe_face.get_function_values(current_solution, values_here);
    fe_face.get_function_gradients(current_solution, gradients_here);

    this->boundary_index = boundary_index;

    const double extent1 = cell->measure() / cell->face(face)->measure();
    this->penalty_factor =
      get_penalty_factor(fe_values_tracer.get_fe().degree, extent1, extent1);

    if (properties_manager.field_is_required(field::levelset))
      {
        Assert(
          levelset_function != nullptr,
          ExcMessage(
            "Levelset function is required for tracer assembly, but the level set function is a nullptr"));

        levelset_function->value_list(face_quadrature_points, sdf_face_values);
      }
  }

  /** @brief Reinitializes the content of the scratch regarding the velocity for internal/boundary faces.
   *  The velocity is inherently assumed to have been solved using a CG scheme.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face The face index associated with the cell
   *
   * @param[in] velocity_solution The present value of the solution.
   * there are any).
   */
  template <typename VectorType>
  void
  reinit_face_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const unsigned int                                   &f,
    const VectorType                                     &velocity_solution,
    const Parameters::ALE<dim>                           &ale,
    std::shared_ptr<Functions::ParsedFunction<dim>>       drift_velocity)
  {
    fe_face_values_fd.reinit(velocity_cell, f);

    // BB note : Array could be pre-allocated
    face_velocity_values.resize(face_quadrature_points.size());

    fe_face_values_fd[velocities].get_function_values(velocity_solution,
                                                      face_velocity_values);


    // Add the drift velocity to the velocity to account for tracer drift flux
    // modeling
    Tensor<1, dim> drift_velocity_tensor;
    Vector<double> drift_velocity_vector(dim);

    for (unsigned int q = 0; q < face_quadrature_points.size(); ++q)
      {
        drift_velocity->vector_value(face_quadrature_points[q],
                                     drift_velocity_vector);
        for (unsigned int d = 0; d < dim; ++d)
          drift_velocity_tensor[d] = drift_velocity_vector[d];

        face_velocity_values[q] += drift_velocity_tensor;
      }


    if (!ale.enabled())
      return;

    // ALE enabled, so extract the ALE velocity and subtract it from the
    // velocity obtained from the fluid dynamics
    Tensor<1, dim>                                  velocity_ale;
    std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
      ale.velocity;
    Vector<double> velocity_ale_vector(dim);

    for (unsigned int q = 0; q < face_quadrature_points.size(); ++q)
      {
        velocity_ale_function->vector_value(face_quadrature_points[q],
                                            velocity_ale_vector);
        for (unsigned int d = 0; d < dim; ++d)
          velocity_ale[d] = velocity_ale_vector[d];

        face_velocity_values[q] -= velocity_ale;
      }
  }

  /** @brief Calculates the physical properties at a face
   */
  void
  calculate_face_physical_properties()
  {
    const auto diffusivity_model = properties_manager.get_tracer_diffusivity();

    // BB note : Array could be pre-allocated
    tracer_diffusivity_face.resize(face_quadrature_points.size());

    // If the trace diffusivity depends on the level set, take this into account
    if (properties_manager.field_is_required(field::levelset))
      set_field_vector(field::levelset, sdf_face_values, fields);


    diffusivity_model->vector_value(fields, tracer_diffusivity_face);
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
  std::vector<double>                  tracer_diffusivity_face;
  std::vector<double>                  tracer_diffusivity_0;
  std::vector<double>                  tracer_diffusivity_1;

  // FEValues for the Tracer problem
  FEValues<dim>          fe_values_tracer;
  FEInterfaceValues<dim> fe_interface_values_tracer;

  unsigned int n_dofs;
  unsigned int n_interface_dofs;
  unsigned int n_q_points;
  double       cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;
  std::vector<Point<dim>> face_quadrature_points;



  // Tracer values
  std::vector<double>              tracer_values;
  std::vector<Tensor<1, dim>>      tracer_gradients;
  std::vector<double>              tracer_laplacians;
  std::vector<std::vector<double>> previous_tracer_values;
  std::vector<Tensor<1, dim>>      previous_tracer_gradients;

  // Tracer values at the faces
  std::vector<double>         values_here;
  std::vector<double>         values_there;
  std::vector<double>         tracer_value_jump;
  std::vector<Tensor<1, dim>> gradients_here;
  std::vector<Tensor<1, dim>> tracer_average_gradient;

  // SIPG (interior faces) or Nitsche (boundary) penalization factor
  // The penalty factor is calculated using a measure of the element size. It is
  // multiplied by the diffusion coefficient to apply the SIPG method
  double penalty_factor;

  // Bool that defines if the selected face is a dirichlet/outlet boundary
  unsigned int boundary_index;

  // Solid signed distance function
  std::vector<double> sdf_values;
  std::vector<double> sdf_face_values;


  // Source term
  std::vector<double> source;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi;
  std::vector<std::vector<double>>         laplacian_phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  FEFaceValues<dim>           fe_face_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<Tensor<1, dim>> face_velocity_values;
};

#endif
