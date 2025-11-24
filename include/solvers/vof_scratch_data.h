// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_scratch_data_h
#define lethe_vof_scratch_data_h

#include <core/time_integration_utilities.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

/**
 * @brief Store the information required by the assembly procedure
 * for a VOF free surface equation. Consequently, this class
 * calculates the phase values (values, gradients, laplacians) and the shape
 * method (values, gradients, laplacians) at all the gauss points for all
 * degrees of freedom and stores it into arrays.
 * This class serves as a separation between the evaluation at the gauss point
 * of the variables of interest and their use in the assembly, which is carried
 * out by the assembler methods.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 **/

template <int dim>
class VOFScratchData : public PhysicsScratchDataBase
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
   * simulation. This is used to extrapolate velocity solutions in time
   * for transient simulation.
   *
   * @param properties_manager Manager to calculate the physical properties.
   *
   * @param fe_vof The FESystem used to solve the VOF equations.
   *
   * @param quadrature The quadrature to use for the assembly within the cells.
   *
   * @param face_quadrature The quadrature to use for the assembly on faces.
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes
   * equations are solved.
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations.
   *
   */
  VOFScratchData(const std::shared_ptr<SimulationControl> &simulation_control,
                 const PhysicalPropertiesManager          &properties_manager,
                 const FiniteElement<dim>                 &fe_vof,
                 const Quadrature<dim>                    &quadrature,
                 const Quadrature<dim - 1>                &face_quadrature,
                 const Mapping<dim>                       &mapping,
                 const FiniteElement<dim>                 &fe_fd)
    : simulation_control(simulation_control)
    , properties_manager(properties_manager)
    , fe_values_vof(mapping,
                    fe_vof,
                    quadrature,
                    update_values | update_gradients |
                      update_quadrature_points | update_hessians |
                      update_JxW_values)
    , fe_interface_values_vof(mapping,
                              fe_vof,
                              face_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values | update_normal_vectors)
    , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
    , fe_face_values_fd(mapping, fe_fd, face_quadrature, update_values)
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
   * @param sd The scratch data to be copied
   */
  VOFScratchData(const VOFScratchData<dim> &sd)
    : simulation_control(sd.simulation_control)
    , properties_manager(sd.properties_manager)
    , fe_values_vof(sd.fe_values_vof.get_mapping(),
                    sd.fe_values_vof.get_fe(),
                    sd.fe_values_vof.get_quadrature(),
                    update_values | update_gradients |
                      update_quadrature_points | update_hessians |
                      update_JxW_values)
    , fe_interface_values_vof(sd.fe_interface_values_vof.get_mapping(),
                              sd.fe_interface_values_vof.get_fe(),
                              sd.fe_interface_values_vof.get_quadrature(),
                              update_values | update_quadrature_points |
                                update_JxW_values | update_normal_vectors)
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values | update_gradients)
    , fe_face_values_fd(sd.fe_face_values_fd.get_mapping(),
                        sd.fe_face_values_fd.get_fe(),
                        sd.fe_face_values_fd.get_quadrature(),
                        update_values)
  {
    allocate();
  }


  /** @brief Allocates the memory for the scratch
   *
   * This method allocates the necessary memory for all members of the scratch
   */
  void
  allocate() override;

  /** @brief Reinitialize the content of the scratch
   *
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the FE which is used to fill the
   * FeValues.
   *
   * @param[in] current_solution The present value of the solution for the VOF.
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions)
  {
    fe_values_vof.reinit(cell);
    this->quadrature_points = fe_values_vof.get_quadrature_points();
    auto &fe_vof            = fe_values_vof.get_fe();

    // Compute cell diameter
    double cell_measure =
      compute_cell_measure_with_JxW(this->fe_values_vof.get_JxW_values());
    this->cell_size = compute_cell_diameter<dim>(cell_measure, fe_vof.degree);

    fe_values_vof.get_function_values(current_solution,
                                      this->present_phase_values);
    fe_values_vof.get_function_gradients(current_solution,
                                         this->phase_gradients);
    fe_values_vof.get_function_laplacians(current_solution,
                                          this->phase_laplacians);
    fe_values_vof.get_function_gradients(previous_solutions[0],
                                         this->previous_phase_gradients);


    // Gather previous vof values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        fe_values_vof.get_function_values(previous_solutions[p],
                                          this->previous_phase_values[p]);
      }

    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->JxW[q] = fe_values_vof.JxW(q);

        for (unsigned int k = 0; k < this->n_dofs; ++k)
          {
            // Shape function
            this->phi[q][k]           = fe_values_vof.shape_value(k, q);
            this->grad_phi[q][k]      = fe_values_vof.shape_grad(k, q);
            this->hess_phi[q][k]      = fe_values_vof.shape_hessian(k, q);
            this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
          }
      }
  }


  /** @brief Reinitialize the content of the scratch for the internal faces. This is only used for the DG assemblers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face_no The face index associated with the cell
   *
   * @param[in] sub_face_no The subface index associated with the face
   *
   * @param[in] neighbor_cell The neighboring cell
   *
   * @param[in] neighbor_face_no The face index associated with the neighboring
   * cell
   *
   * @param[in] neighbor_sub_face_no The subface index associated with the
   * neighboring cell
   *
   * @param[in] current_solution The present value of the solution.
   */
  template <typename VectorType>
  void
  reinit_internal_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                   &face_no,
    const unsigned int                                   &sub_face_no,
    const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
    const unsigned int                                   &neighbor_face_no,
    const unsigned int                                   &neighbor_sub_face_no,
    const VectorType                                     &current_solution)
  {
    fe_interface_values_vof.reinit(cell,
                                   face_no,
                                   sub_face_no,
                                   neighbor_cell,
                                   neighbor_face_no,
                                   neighbor_sub_face_no);
    face_quadrature_points = fe_interface_values_vof.get_quadrature_points();

    n_interface_dofs = fe_interface_values_vof.n_current_interface_dofs();

    // BB TODO : Preallocate memory here
    values_here.resize(face_quadrature_points.size());
    values_there.resize(face_quadrature_points.size());
    phase_value_jump.resize(face_quadrature_points.size());

    fe_interface_values_vof.get_fe_face_values(0).get_function_values(
      current_solution, values_here);
    fe_interface_values_vof.get_fe_face_values(1).get_function_values(
      current_solution, values_there);

    fe_interface_values_vof.get_jump_in_function_values(current_solution,
                                                        phase_value_jump);
  }

  /** @brief Reinitialize the velocity, calculated by the Fluid Dynamics
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell for which the velocity is reinitialized
   * This cell must be compatible with the Fluid Dynamics FE
   *
   * @param current_solution The present value of the solution for \f$[u,p]\f$
   *
   * @param previous_solutions Vector of \f$n\f$ @p VectorType containers of
   * previous fluid dynamic solutions (\f$[u,p]\f$). \f$n\f$ depends on the BDF
   * scheme selected for time-stepping.
   *
   * @param ale ALE parameters
   *
   */

  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType              &current_solution,
                  const std::vector<VectorType> &previous_solutions,
                  const Parameters::ALE<dim>    &ale)
  {
    fe_values_fd.reinit(cell);

    fe_values_fd[velocities_fd].get_function_values(current_solution,
                                                    velocity_values);
    fe_values_fd[velocities_fd].get_function_gradients(
      current_solution, velocity_gradient_values);

    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        this->velocity_divergences[q] =
          trace(this->velocity_gradient_values[q]);
      }

    // Gather previous velocity values
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        fe_values_fd[velocities_fd].get_function_values(
          previous_solutions[p], this->previous_velocity_values[p]);
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
            for (int d = 0; d < dim; ++d)
              velocity_ale[d] = velocity_ale_vector[d];

            velocity_values[q] -= velocity_ale;

            for (unsigned int p = 0; p < previous_solutions.size(); ++p)
              {
                this->previous_velocity_values[p][q] -= velocity_ale;
              }
          }
      }

    // Extrapolate velocity to t+dt using the BDF scheme if the simulation is
    // transient
    const auto method = this->simulation_control->get_assembly_method();
    if (time_stepping_is_bdf(method))
      {
        // Extrapolate velocity
        std::vector<double> time_vector =
          this->simulation_control->get_simulation_times();
        bdf_extrapolate(time_vector,
                        this->previous_velocity_values,
                        number_of_previous_solutions(method),
                        this->velocity_values);
      }
  }


  /** @brief Reinitialize the content of the scratch regarding the velocity for internal/boundary faces.
   *  The velocity is inherently assumed to have been solved using a CG scheme.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face_no The face index associated with the cell
   *
   * @param[in] velocity_solution The present value of the velocity solution.
   */
  template <typename VectorType>
  void
  reinit_face_velocity(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const unsigned int                                   &face_no,
    const VectorType                                     &velocity_solution,
    const Parameters::ALE<dim>                           &ale)
  {
    fe_face_values_fd.reinit(velocity_cell, face_no);

    // BB note : Array could be pre-allocated
    face_velocity_values.resize(face_quadrature_points.size());

    fe_face_values_fd[velocities_fd].get_function_values(velocity_solution,
                                                         face_velocity_values);

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
        for (int d = 0; d < dim; ++d)
          velocity_ale[d] = velocity_ale_vector[d];

        face_velocity_values[q] -= velocity_ale;
      }
  }

  // For velocity solution extrapolation
  const std::shared_ptr<SimulationControl> simulation_control;

  // Physical properties
  const PhysicalPropertiesManager      properties_manager;
  std::map<field, std::vector<double>> fields;

  // FEValues for the VOF problem
  FEValues<dim>          fe_values_vof;
  FEInterfaceValues<dim> fe_interface_values_vof;
  unsigned int           n_dofs;
  unsigned int           n_interface_dofs;
  unsigned int           n_q_points;
  double                 cell_size;

  // Quadrature
  Table<1, double>        JxW;
  std::vector<Point<dim>> quadrature_points;
  std::vector<Point<dim>> face_quadrature_points;

  // VOF values
  std::vector<double>         present_phase_values;
  std::vector<Tensor<1, dim>> phase_gradients;
  std::vector<Tensor<1, dim>> previous_phase_gradients;

  std::vector<double>              phase_laplacians;
  std::vector<std::vector<double>> previous_phase_values;

  // VOF values at the faces
  std::vector<double> values_here;
  std::vector<double> values_there;
  std::vector<double> phase_value_jump;

  // Shape functions
  Table<2, double>         phi;
  Table<2, Tensor<1, dim>> grad_phi;
  Table<2, Tensor<2, dim>> hess_phi;
  Table<2, double>         laplacian_phi;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValues<dim>     fe_values_fd;
  FEFaceValues<dim> fe_face_values_fd;

  FEValuesExtractors::Vector velocities_fd;
  // This FEValues must be instantiated for the velocity
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<Tensor<2, dim>>              velocity_gradient_values;
  std::vector<double>                      velocity_divergences;

  // Face velocity value for DG
  std::vector<Tensor<1, dim>> face_velocity_values;
};

#endif
