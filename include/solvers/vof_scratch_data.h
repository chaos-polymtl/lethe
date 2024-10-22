// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_vof_scratch_data_h
#define lethe_vof_scratch_data_h

#include <core/multiphysics.h>
#include <core/time_integration_utilities.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>

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
   * @param fe_vof The FESystem used to solve the VOF equations
   *
   * @param quadrature The quadrature to use for the assembly
   *
   * @param mapping The mapping of the domain in which the Navier-Stokes
   * equations are solved
   *
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  VOFScratchData(const std::shared_ptr<SimulationControl> &simulation_control,
                 const PhysicalPropertiesManager          &properties_manager,
                 const FiniteElement<dim>                 &fe_vof,
                 const Quadrature<dim>                    &quadrature,
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
    , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
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
    , fe_values_fd(sd.fe_values_fd.get_mapping(),
                   sd.fe_values_fd.get_fe(),
                   sd.fe_values_fd.get_quadrature(),
                   update_values | update_gradients)
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
            for (unsigned int d = 0; d < dim; ++d)
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
    if (is_bdf(method))
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

  // For velocity solution extrapolation
  const std::shared_ptr<SimulationControl> simulation_control;

  // Physical properties
  const PhysicalPropertiesManager      properties_manager;
  std::map<field, std::vector<double>> fields;

  // FEValues for the VOF problem
  FEValues<dim> fe_values_vof;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  std::vector<double>     JxW;
  std::vector<Point<dim>> quadrature_points;

  // VOF values
  std::vector<double>         present_phase_values;
  std::vector<Tensor<1, dim>> phase_gradients;
  std::vector<Tensor<1, dim>> previous_phase_gradients;

  std::vector<double>              phase_laplacians;
  std::vector<std::vector<double>> previous_phase_values;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;
  std::vector<std::vector<Tensor<2, dim>>> hess_phi;
  std::vector<std::vector<double>>         laplacian_phi;


  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValues<dim> fe_values_fd;

  FEValuesExtractors::Vector velocities_fd;
  // This FEValues must be instantiated for the velocity
  std::vector<Tensor<1, dim>>              velocity_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_velocity_values;
  std::vector<Tensor<2, dim>>              velocity_gradient_values;
  std::vector<double>                      velocity_divergences;
};

/**
 * @brief Store the information required by the assembly procedure
 * for a VOF phase gradient projection in the L2 space.
 * Computes the phase gradient values and gradients and the shape
 * function values and gradients at all the quadrature points for all
 * degrees of freedom of a given cell and stores it into arrays.
 * Serves as a link between the evaluation at the quadrature point
 * of the variables of interest and their use in the assembly.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup solvers
 **/
template <int dim>
class VOFPhaseGradientProjectionScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Constructor of the scratch data object.
   * Creates the fe_values that will be used to fill the member variables.
   * Allocates the memory necessary for all member variables, but does not
   * compute anything since that must be done at the cell level.
   *
   * @param[in] fe_phase_gradient_projection FiniteElement object used for
   * phase fraction gradient L2 projection.
   *
   * @param[in] quadrature Quadrature rule used for the assembly of the matrix
   * and the right-hand side.
   *
   * @param[in] mapping Mapping of the domain used when solving the VOF
   * equation.
   *
   * @param[in] fe_vof FiniteElement object used for VOF.
   */
  VOFPhaseGradientProjectionScratchData(
    const FiniteElement<dim> &fe_phase_gradient_projection,
    const Quadrature<dim>    &quadrature,
    const Mapping<dim>       &mapping,
    const FiniteElement<dim> &fe_vof)
    : fe_values_phase_gradient_projection(mapping,
                                          fe_phase_gradient_projection,
                                          quadrature,
                                          update_values | update_gradients |
                                            update_JxW_values)
    , fe_values_vof(mapping, fe_vof, quadrature, update_gradients)
  {
    allocate();
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * The constructor only uses the other scratch to build the FEValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param[in] sd The scratch data to be copied.
   */
  VOFPhaseGradientProjectionScratchData(
    const VOFPhaseGradientProjectionScratchData<dim> &sd)
    : fe_values_phase_gradient_projection(
        sd.fe_values_phase_gradient_projection.get_mapping(),
        sd.fe_values_phase_gradient_projection.get_fe(),
        sd.fe_values_phase_gradient_projection.get_quadrature(),
        update_values | update_gradients | update_JxW_values)
    , fe_values_vof(sd.fe_values_vof.get_mapping(),
                    sd.fe_values_vof.get_fe(),
                    sd.fe_values_vof.get_quadrature(),
                    update_gradients)
  {
    allocate();
  }

  /**
   * @brief Default destructor.
   */
  ~VOFPhaseGradientProjectionScratchData() = default;

  /**
   * @brief Allocate the memory necessary memory for all members of the
   * scratch.
   */
  void
  allocate() override;

  /**
   * @brief Reinitialize the content of the object for a given cell.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the reinitialization is carried out.
   *
   * @param[in] current_solution Present values of the solution for the
   * projected phase fraction gradient.
   */
  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution)
  {
    // Reinitialize FEValues
    this->fe_values_phase_gradient_projection.reinit(cell);
    auto &fe_phase_gradient_projection =
      this->fe_values_phase_gradient_projection.get_fe();

    // Get cell measure
    this->JxW = this->fe_values_phase_gradient_projection.get_JxW_values();
    double cell_measure = compute_cell_measure_with_JxW(this->JxW);
    this->cell_size =
      compute_cell_diameter<dim>(cell_measure,
                                 fe_phase_gradient_projection.degree);

    // Gather present solutions
    this->fe_values_phase_gradient_projection[this->phase_fraction_gradients]
      .get_function_values(current_solution,
                           this->present_phase_gradient_projection_values);
    this->fe_values_phase_gradient_projection[this->phase_fraction_gradients]
      .get_function_gradients(
        current_solution, this->present_phase_gradient_projection_gradients);

    // Shape functions
    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        for (unsigned int k = 0; k < this->n_dofs; ++k)
          {
            this->phi[q][k] = this
                                ->fe_values_phase_gradient_projection
                                  [this->phase_fraction_gradients]
                                .value(k, q);
            this->grad_phi[q][k] = this
                                     ->fe_values_phase_gradient_projection
                                       [this->phase_fraction_gradients]
                                     .gradient(k, q);
          }
      }
  }

  /**
   * @brief Reinitialize the cell with the VOF solution calculated with the VOF
   * advection equation.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the assembly is carried out.
   *
   * @param[in] current_solution Present value of the solution for the
   * VOF phase fraction.
   */
  template <typename VectorType>
  void
  reinit_vof(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType &current_solution)
  {
    // Reinitialize FEValues to cell
    this->fe_values_vof.reinit(cell);

    // Gather present VOF phase fraction gradients
    this->fe_values_vof.get_function_gradients(
      current_solution, present_filtered_vof_phase_gradients);
  }

  // FEValues for the VOF phase gradient projection problem
  FEValues<dim>              fe_values_phase_gradient_projection;
  unsigned int               n_dofs;
  unsigned int               n_q_points;
  double                     cell_size;
  FEValuesExtractors::Vector phase_fraction_gradients;

  // Jacobi determinant times the quadrature weights
  std::vector<double> JxW;

  // Phase gradient projection values and gradients
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values;
  std::vector<Tensor<2, dim>> present_phase_gradient_projection_gradients;

  // VOF phase gradient values
  std::vector<Tensor<1, dim>> present_filtered_vof_phase_gradients;

  // Shape functions
  std::vector<std::vector<Tensor<1, dim>>> phi;
  std::vector<std::vector<Tensor<2, dim>>> grad_phi;

  // FEValues for the VOF problem
  FEValues<dim> fe_values_vof;
};


/**
 * @brief Store the information required by the assembly procedure
 * for a VOF curvature projection in the L2 space.
 * Computes the curvature values and gradients and the shape
 * function values and gradients at all the quadrature points for all
 * degrees of freedom of a given cell and stores it into arrays.
 * Serves as a link between the evaluation at the quadrature point
 * of the variables of interest and their use in the assembly.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @ingroup solvers
 **/
template <int dim>
class VOFCurvatureProjectionScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Constructor of the scratch data object.
   * Creates the fe_values that will be used to fill the member variables.
   * Allocates the memory necessary for all member variables, but does not
   * compute anything since that must be done at the cell level.
   *
   * @param[in] fe_curvature_projection FiniteElement object used for
   * curvature L2 projection.
   *
   * @param[in] quadrature Quadrature rule used for the assembly of the matrix
   * and the right-hand side.
   *
   * @param[in] mapping Mapping of the domain used when solving the VOF
   * equation.
   *
   * @param[in] fe_phase_gradient_projection FiniteElement object used for
   * phase fraction gradient L2 projection.
   */
  VOFCurvatureProjectionScratchData(
    const FiniteElement<dim> &fe_curvature_projection,
    const Quadrature<dim>    &quadrature,
    const Mapping<dim>       &mapping,
    const FiniteElement<dim> &fe_phase_gradient_projection)
    : fe_values_curvature_projection(mapping,
                                     fe_curvature_projection,
                                     quadrature,
                                     update_values | update_gradients |
                                       update_JxW_values)
    , fe_values_phase_gradient_projection(mapping,
                                          fe_phase_gradient_projection,
                                          quadrature,
                                          update_values)
  {
    allocate();
  }

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * The constructor only uses the other scratch to build the FEValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param[in] sd The scratch data to be copied.
   */
  VOFCurvatureProjectionScratchData(
    const VOFCurvatureProjectionScratchData<dim> &sd)
    : fe_values_curvature_projection(
        sd.fe_values_curvature_projection.get_mapping(),
        sd.fe_values_curvature_projection.get_fe(),
        sd.fe_values_curvature_projection.get_quadrature(),
        update_values | update_gradients | update_JxW_values)
    , fe_values_phase_gradient_projection(
        sd.fe_values_phase_gradient_projection.get_mapping(),
        sd.fe_values_phase_gradient_projection.get_fe(),
        sd.fe_values_phase_gradient_projection.get_quadrature(),
        update_values)
  {
    allocate();
  }

  /**
   * @brief Default destructor.
   */
  ~VOFCurvatureProjectionScratchData() = default;
  /**
   * @brief Allocate the memory necessary memory for all members of the
   * scratch.
   */
  void
  allocate() override;

  /**
   * @brief Reinitialize the content of the object for a given cell.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the reinitialization is carried out.
   *
   * @param[in] current_solution Present values of the solution for the
   * curvature projection.
   */
  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution)
  {
    // Reinitialize FEValues
    this->fe_values_curvature_projection.reinit(cell);
    auto &fe_curvature_projection =
      this->fe_values_curvature_projection.get_fe();
    // Get cell measure
    this->JxW           = this->fe_values_curvature_projection.get_JxW_values();
    double cell_measure = compute_cell_measure_with_JxW(this->JxW);
    this->cell_size =
      compute_cell_diameter<dim>(cell_measure, fe_curvature_projection.degree);
    // Gather present solutions
    this->fe_values_curvature_projection.get_function_values(
      current_solution, this->present_curvature_projection_values);
    this->fe_values_curvature_projection.get_function_gradients(
      current_solution, this->present_curvature_projection_gradients);
    // Shape functions
    for (unsigned int q = 0; q < this->n_q_points; ++q)
      {
        for (unsigned int k = 0; k < this->n_dofs; ++k)
          {
            this->phi[q][k] =
              this->fe_values_curvature_projection.shape_value(k, q);
            this->grad_phi[q][k] =
              this->fe_values_curvature_projection.shape_grad(k, q);
          }
      }
  }

  /**
   * @brief Reinitialize the cell with the VOF projected phase fraction gradient
   * solution.
   *
   * @tparam VectorType Vector type used for the solvers.
   *
   * @param[in] cell Cell over which the assembly is carried out.
   *
   * @param[in] current_solution Present values of the VOF projected phase
   * fraction gradient solution.
   */
  template <typename VectorType>
  void
  reinit_projected_phase_gradient(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution)
  {
    // Reinitialize FEValues to cell
    this->fe_values_phase_gradient_projection.reinit(cell);
    // Gather present VOF phase fraction gradients
    this->fe_values_phase_gradient_projection[this->phase_fraction_gradients]
      .get_function_values(current_solution,
                           present_phase_gradient_projection_values);
  }

  // FEValues for the VOF curvature projection problem
  FEValues<dim> fe_values_curvature_projection;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Jacobi determinant times the quadrature weights
  std::vector<double> JxW;

  // Curvature projection values and gradients
  std::vector<double>         present_curvature_projection_values;
  std::vector<Tensor<1, dim>> present_curvature_projection_gradients;

  // Projected phase gradient values
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values;
  FEValuesExtractors::Vector  phase_fraction_gradients;

  // Shape functions
  std::vector<std::vector<double>>         phi;
  std::vector<std::vector<Tensor<1, dim>>> grad_phi;

  // FEValues for the VOF phase gradient projection problem
  FEValues<dim> fe_values_phase_gradient_projection;
};

#endif
