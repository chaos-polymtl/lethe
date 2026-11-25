// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_heat_transfer_scratch_data_h
#define lethe_heat_transfer_scratch_data_h

#include <core/physical_property_model.h>
#include <core/vector.h>

#include <solvers/cls_filter.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/fe/fe_system.h>

using namespace dealii;

/**
 * @brief Class that stores the information required by the assembly procedure
 * for a heat transfer advection-diffusion equation. Consequently, this class
 * calculates the heat transfer (values, gradients, laplacians) and the shape
 * function (values, gradients, laplacians) at all the gauss points for all
 * degrees of freedom and stores it into arrays. Additionally, the user can
 * request that this class gathers additional fields for physics which are
 * coupled to the heat transfer equation, such as the velocity which is
 * required. This class serves as a separation between the evaluation at the
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
class HeatTransferScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evaluation, since this needs to be done at the cell level.
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
  HeatTransferScratchData(const PhysicalPropertiesManager &properties_manager,
                          const FiniteElement<dim>        &fe_ht,
                          const Quadrature<dim>           &quadrature,
                          const Mapping<dim>              &mapping,
                          const FiniteElement<dim>        &fe_fd,
                          const Quadrature<dim - 1>       &face_quadrature,
                          const double                     delta_T_ref);

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism, it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd The scratch data
   */
  HeatTransferScratchData(const HeatTransferScratchData<dim> &sd);

  /** @brief Allocates the memory for the scratch
   *
   * This method allocates the necessary memory for all members of the scratch
   *
   */
  void
  allocate() override;

  /** @brief Reinitializes the content of the scratch.
   *
   * Using the FeValues and the content of the solutions, previous solutions,
   * fills all of the class member of the scratch.
   *
   * @tparam VectorType The Vector type used for the solvers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the FE which is used to fill the
   * FeValues.
   *
   * @param[in] current_solution The present value of the solution for the Heat
   * Transfer.
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   * @param[in] source_function The function describing the Heat Transfer source
   * term.
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<GlobalVectorType> &previous_solutions,
         Function<dim>                       *source_function);

  /**
   * @brief Reinitialize cell's velocity and pressure values with values
   * calculated in the Fluid Dynamics (FD) physic while also taking into account
   * ALE for the velocity.
   *
   * @tparam VectorType Type of vector used for storing fluid dynamics
   * solutions.
   *
   * @param[in] cell Pointer to an active cell for which the velocity and
   * pressure are reinitialized.This cell must be compatible with the FD
   * FiniteElement.
   *
   * @param[in] current_solution Solution vector from FD containing velocity
   * components and pressure values.
   *
   * @param[in] ale The ALE parameters which include the ALE function
   *
   */

  template <typename VectorType>
  void
  reinit_fluid_dynamics(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution,
    const Parameters::ALE<dim>                           &ale);

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
   * @brief enable_cls Enables the collection of the CLS data by the scratch.
   *
   * @param fe FiniteElement associated with the CLS.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly.
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param phase_filter_parameters Parameters for phase indicator filtering.
   */

  void
  enable_cls(const FiniteElement<dim>          &fe,
             const Quadrature<dim>             &quadrature,
             const Mapping<dim>                &mapping,
             const Parameters::CLS_PhaseFilter &phase_filter_parameters);

  /**
   * @brief enable_cls Enables the collection of the CLS data by the scratch - function overload used in the copy constructor of HeatTransferScratchData
   *
   * @param fe FiniteElement associated with the CLS.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param filter Filter that is applied on the phase indicator
   */

  void
  enable_cls(const FiniteElement<dim>                              &fe,
             const Quadrature<dim>                                 &quadrature,
             const Mapping<dim>                                    &mapping,
             const std::shared_ptr<ConservativeLevelSetFilterBase> &filter);

  /** @brief Reinitialize the content of the scratch for CLS.
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the CLS FE and not the
   * Fluid Dynamics FE
   *
   * @param current_solution The present solution for the phase value
   *
   * NB: the previous_solutions are not used in heat_transfer_assemblers,
   * contrary to navier_stokes_assemblers
   */

  template <typename VectorType>
  void
  reinit_cls(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType &current_solution,
             const VectorType &current_filtered_solution);

  /**
   * @brief enable_time_harmonic_maxwell Enables the collection of the
   * time-harmonic Maxwell data by the scratch - function overload used in the
   * copy constructor of HeatTransferScratchData
   *
   * @param fe FiniteElement associated with the time-harmonic Maxwell problem.
   *
   * @param quadrature Quadrature rule of the time-harmonic Maxwell problem
   * assembly
   *
   * @param mapping Mapping used for the time-harmonic Maxwell problem assembly
   */

  void
  enable_time_harmonic_maxwell(const FiniteElement<dim> &fe,
                               const Quadrature<dim>    &quadrature,
                               const Mapping<dim>       &mapping);

  /** @brief Reinitialize the content of the scratch for the time-harmonic
   * Maxwell problem.
   *
   * @tparam VectorType The Vector type used for the solvers
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the time-harmonic Maxwell FE and not the
   * Fluid Dynamics FE
   *
   * @param current_solution The present solution for the electromagnetic field
   * values
   */

  template <typename VectorType>
  void
  reinit_time_harmonic_maxwell(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const VectorType                                     &current_solution);

  /** @brief Calculates the physical properties. This method calculates the
   * physical properties that may be required by the heat transfer problem.
   * Namely the density, specific heat, thermal conductivity and viscosity (for
   * viscous dissipation).
   *
   */
  void
  calculate_physical_properties();


  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  dealii::types::material_id           material_id;
  std::vector<double>                  specific_heat;
  std::vector<double>                  thermal_conductivity;
  std::vector<double>                  density;
  std::vector<double>                  dynamic_viscosity;
  // Gradient of the specific heat with respect to the temperature
  // This is calculated by deriving the specific heat by the temperature
  // (dCp/dT)
  std::vector<double> grad_specific_heat_temperature;

  // Auxiliary property vector for CLS simulations
  std::vector<double> specific_heat_0;
  std::vector<double> thermal_conductivity_0;
  std::vector<double> density_0;
  std::vector<double> dynamic_viscosity_0;
  std::vector<double> grad_specific_heat_temperature_0;

  std::vector<double> specific_heat_1;
  std::vector<double> thermal_conductivity_1;
  std::vector<double> density_1;
  std::vector<double> dynamic_viscosity_1;
  std::vector<double> grad_specific_heat_temperature_1;

  // Physical properties for time-harmonic Maxwell's equations. Note that only
  // the imaginary part of the electric permittivity and magnetic permeability
  // are responsible for the heat generation, so only these are stored here (and
  // conductivity is always real).
  std::vector<double> electric_conductivity;
  std::vector<double> magnetic_permeability_imag;
  std::vector<double> electric_permittivity_imag;

  // FEValues for the HT problem
  FEValues<dim> fe_values_T;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;


  // Quadrature
  Table<1, double>        JxW;
  std::vector<Point<dim>> quadrature_points;

  // Temperature values
  std::vector<double>                      present_temperature_values;
  std::vector<Tensor<1, dim>>              temperature_gradients;
  std::vector<double>                      present_temperature_laplacians;
  std::vector<double>                      present_face_temperature_values;
  std::vector<std::vector<double>>         previous_temperature_values;
  std::vector<std::vector<Tensor<1, dim>>> previous_temperature_gradients;
  double                                   global_delta_T_ref;

  // Shape functions and gradients
  Table<2, double>         phi_T;
  Table<2, Tensor<1, dim>> grad_phi_T;
  Table<2, Tensor<2, dim>> hess_phi_T;
  Table<2, double>         laplacian_phi_T;

  // Source term
  std::vector<double> source;

  /**
   * Scratch component for the CLS auxiliary physics
   */
  bool                        gather_cls;
  unsigned int                n_dofs_cls;
  std::vector<double>         filtered_phase_values;
  std::vector<Tensor<1, dim>> filtered_phase_gradient_values;
  std::vector<Tensor<1, dim>> phase_gradient_values;

  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_cls;
  std::shared_ptr<ConservativeLevelSetFilterBase>
    filter; // Phase indicator filter

  /**
   * Scratch component for the time-harmonic Maxwell auxiliary physics
   */
  bool                        gather_time_harmonic_maxwell;
  FEValuesExtractors::Vector  extractor_E_real;
  FEValuesExtractors::Vector  extractor_E_imag;
  FEValuesExtractors::Vector  extractor_H_real;
  FEValuesExtractors::Vector  extractor_H_imag;
  std::vector<Tensor<1, dim>> electric_field_real_values;
  std::vector<Tensor<1, dim>> electric_field_imag_values;
  std::vector<Tensor<1, dim>> magnetic_field_real_values;
  std::vector<Tensor<1, dim>> magnetic_field_imag_values;

  // Stored as shared_ptr because only instantiated when EM is active
  std::shared_ptr<FEValues<dim>> fe_values_time_harmonic_maxwell;

  /**
   * Scratch component for the Navier-Stokes components
   */
  const FEValuesExtractors::Vector velocities;
  const FEValuesExtractors::Scalar pressure;
  // This FEValues must be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<Tensor<2, dim>> velocity_gradient_values;
  std::vector<double>         shear_rate_values;
  /// Pressure field for physical property models requiring them
  std::vector<double> pressure_values;

  // Scratch for the face boundary condition
  FEFaceValues<dim> fe_face_values_ht;
  Table<2, double>  face_JxW;

  unsigned int n_faces;
  unsigned int n_faces_q_points;

  // Is boundary cell indicator
  bool                      is_boundary_cell;
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;


  // First is face number, second quadrature point, third DOF
  Table<3, double> phi_face_T;
  // First vector is face number, second quadrature point
  std::vector<std::vector<double>> temperature_face_value;
};

#endif
