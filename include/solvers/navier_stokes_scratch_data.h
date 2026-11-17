// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_navier_stokes_scratch_data_h
#define lethe_navier_stokes_scratch_data_h

#include "core/parameters_cfd_dem.h"
#include <core/ale.h>
#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/lethe_grid_tools.h>
#include <core/parameters.h>
#include <core/physical_property_model.h>
#include <core/sdirk_stage_data.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_filter.h>
#include <solvers/cls_filter.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_system.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Class that stores the information required by the assembly procedure
 * for a Navier-Stokes equation. Consequently, this class calculates
 * the velocity (values, gradients, laplacians) and the shape function
 * (values, gradients, laplacians) at all the gauss points for all degrees
 * of freedom and stores it into arrays. Additionally, the use can request
 * that this class gathers additional fields for physics which are coupled
 * to the Navier-Stokes equation, such as the CLS. This class
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
  // Store the curl using deal.II's own view-dependent result type. This keeps
  // the scratch data compatible with the deal.II versions used in CI, where
  // the 2D curl may be exposed either as a scalar or as Tensor<1,1>.
  using velocity_curl_type_2d =
    typename FEValuesViews::Vector<2>::template solution_curl_type<double>;
  using velocity_curl_type_3d =
    typename FEValuesViews::Vector<3>::template solution_curl_type<double>;

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
   * @param face_quadrature The quadrature to use on faces for the assembly
   */
  NavierStokesScratchData(
    const std::shared_ptr<SimulationControl> &simulation_control,
    const PhysicalPropertiesManager          &properties_manager,
    const FESystem<dim>                      &fe,
    const Quadrature<dim>                    &quadrature,
    const Mapping<dim>                       &mapping,
    const Quadrature<dim - 1>                &face_quadrature);

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd The scratch data to be copied
   */
  NavierStokesScratchData(const NavierStokesScratchData<dim> &sd);

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
         const std::vector<VectorType>        &previous_solutions,
         const VectorType                     &sum_over_previous_stages,
         const std::shared_ptr<Function<dim>> &forcing_function,
         const Tensor<1, dim>                 &beta_force,
         const double                          pressure_scaling_factor);

  /**
   * @brief enable_cls Enables the collection of the CLS data by the scratch
   *
   * @param fe FiniteElement associated with the CLS.
   *
   * @param quadrature Quadrature rule of the Navier-Stokes problem assembly
   *
   * @param mapping Mapping used for the Navier-Stokes problem assembly
   *
   * @param phase_filter_parameters Parameters for phase indicator filtering
   */

  void
  enable_cls(const FiniteElement<dim>          &fe,
             const Quadrature<dim>             &quadrature,
             const Mapping<dim>                &mapping,
             const Parameters::CLS_PhaseFilter &phase_filter_parameters);

  /**
   * @brief enable_cls Enables the collection of the CLS data by the scratch - function overload used in the copy constructor of NavierStokesScratchData
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

  void
  enable_projected_phase_indicator_gradient(
    const FiniteElement<dim> &fe_projected_phase_indicator_gradient,
    const Quadrature<dim>    &quadrature,
    const Mapping<dim>       &mapping);

  void
  enable_curvature(const FiniteElement<dim> &fe_curvature,
                   const Quadrature<dim>    &quadrature,
                   const Mapping<dim>       &mapping);

  /**
   * @brief Reinitialize the content of the scratch for the cls
   *
   * @param cell The cell over which the assembly is being carried.
   * This cell must be compatible with the CLS FE and not the
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
  reinit_cls(const typename DoFHandler<dim>::active_cell_iterator &cell,
             const VectorType              &current_solution,
             const VectorType              &current_filtered_solution,
             const std::vector<VectorType> &previous_solutions);

  template <typename VectorType>
  void
  reinit_projected_phase_indicator_gradient(
    const typename DoFHandler<dim>::active_cell_iterator
                     &projected_phase_indicator_gradient_cell,
    const VectorType &current_projected_phase_indicator_gradient_solution);

  template <typename VectorType>
  void
  reinit_curvature(
    const typename DoFHandler<dim>::active_cell_iterator &curvature_cell,
    const VectorType &current_curvature_solution);

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
    const std::vector<VectorType>                        &previous_solutions);

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
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <typename PropertiesIndex>
  void
  reinit_particle_fluid_forces();

  /**
   * @brief Extracts the velocity of the particles and calculates their total volume
   * in the cell. Both outcomes are stored as members of the scratch data class.
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  template <typename PropertiesIndex>
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
    const Parameters::DragCoupling &drag_coupling);

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
    const VectorType &void_fraction_solution);

  /**
   * @brief Calculates the properties of the fluid at the locations of the particles.
   * At the moment, only constant properties within the same fluid are
   * supported. When two fluids are present and CLS is used, the properties are
   * calculated based on the filtered CLS solution interpolated at the
   * location of the particles. These properties are used in the forces
   * calculations in the VANS equations.
   */

  void
  calculate_fluid_properties_at_particle_location();

  /**
   * @brief Calculates the velocity of the fluid relative to that of the particle and
   * the particle Reynolds number at the location of the particles. These will
   * be used in forces calculations in the vans equations.
   *
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */

  template <typename PropertiesIndex>
  void
  calculate_force_parameters_at_particle_location();

  /**
   * @brief Interpolates the filtered CLS solution at the location of the particles.
   * The latter values are used in calculating the density and viscosity of the
   * fluid at the particles' locations when CLS is used.
   *
   * @param[in] q_particles_location Quadrature type object that contains the
   * location of the particles relative to the cell's frame of reference.
   *
   * @param[in] phase_cell The active cell associated with the CLS DoFHandler
   *
   * @param[in] current_filtered_solution The present value of the filtered CLS
   * solution
   */

  template <typename VectorType>
  void
  calculate_cls_at_particle_location(
    const Quadrature<dim>                                &q_particles_location,
    const typename DoFHandler<dim>::active_cell_iterator &phase_cell,
    const VectorType &current_filtered_solution);

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

  template <typename PropertiesIndex, typename VectorType>
  void
  reinit_particle_fluid_interactions(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
    const VectorType                      &present_velocity_pressure_solution,
    const VectorType                      &previous_velocity_pressure_solution,
    const VectorType                      &void_fraction_solution,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::DragCoupling        &drag_coupling);

  /**
   * @brief Calculates the variables needed to compute the particle fluid interactions
   * in the VANS equations.This version of the function is used when CLS is used
   * with CFD-DEM.
   *
   * @param[in] velocity_cell The active cell associated with the velocity and
   * pressure DoFHandler.
   *
   * @param[in] void_fraction_cell The active cell associated with the void
   * fraction DoFHandler.
   *
   * @param[in] phase_cell The active cell associated with the CLS DoFHandler.
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
   * @param[in] current_filtered_CLS_solution The present value of the CLS
   * solution.
   *
   * @param[in] drag_coupling Indicator for the type of coupling that is
   * used. This parameter essentially decides which velocity_pressure solution
   * is interpolated at the location of the particles.
   */
  template <typename PropertiesIndex, typename VectorType>
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
    const VectorType                      &current_filtered_CLS_solution);

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
    const std::vector<VectorType>                        &previous_solutions);

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
    const VectorType &current_filtered_solution);

  /**
   * @brief enable_ale Enables the calculation of the ALE velocity
   */
  void
  enable_ale();

  /**
   * @brief Renitialize the content of the scratch data for ALE
   *
   * @param[in] ale ALE parameters
   */
  void
  reinit_ale(const Parameters::ALE<dim> &ale);

  /**
   * @brief enable_mortar Enables the calculation of ALE velocity within the mortar feature
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
   */
  void
  reinit_mortar(const typename DoFHandler<dim>::active_cell_iterator &cell,
                const Parameters::Mortar<dim> &mortar_parameters);

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
    const Parameters::DragCoupling &drag_coupling);

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
  std::vector<double>              density_thermal_expansion;
  std::vector<double>              grad_kinematic_viscosity_shear_rate;
  std::vector<std::vector<double>> previous_density;

  // Pressure scaling factor to facilitate different scales between velocity and
  // pressure
  double pressure_scaling_factor;

  // For CLS and CH simulations. Present properties for fluid 0 and 1.
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
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<double>         velocity_divergences;
  std::vector<Tensor<2, dim>> velocity_gradients;
  std::vector<Tensor<1, dim>> velocity_laplacians;
  std::vector<Tensor<3, dim>> velocity_hessians;
  std::vector<Tensor<1, dim>> velocity_gradient_divergence;
  // The gradient of the velocity divergence
  std::vector<Tensor<1, dim>>              advective_velocity;
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
  Table<2, Tensor<1, dim>> gradient_divergence_phi_u;
  Table<2, Tensor<2, dim>> grad_phi_u;
  Table<2, double>         phi_p;
  Table<2, Tensor<1, dim>> grad_phi_p;

  /**
   * Scratch component for the CLS auxiliary physics
   */
  bool                             gather_cls;
  unsigned int                     n_dofs_cls;
  std::vector<double>              phase_values;
  std::vector<double>              filtered_phase_values;
  std::vector<double>              filtered_phase_values_at_particle_location;
  std::vector<std::vector<double>> previous_phase_values;
  std::vector<Tensor<1, dim>>      filtered_phase_gradient_values;
  std::vector<Tensor<1, dim>>      phase_gradient_values;
  // This is stored as a shared_ptr because it is only instantiated when needed
  std::shared_ptr<FEValues<dim>> fe_values_cls;
  // Used to evaluate the filtered CLS solution at particle locations
  // (arbitrary points within the cell); stored as a shared_ptr because it is
  // only instantiated when needed
  std::shared_ptr<FEPointEvaluation<1, dim>> point_evaluator_cls;
  std::shared_ptr<ConservativeLevelSetFilterBase>
    filter; // Phase indicator filter

  bool                           gather_projected_phase_indicator_gradient;
  bool                           gather_curvature;
  std::shared_ptr<FEValues<dim>> fe_values_projected_phase_indicator_gradient;
  std::shared_ptr<FEValues<dim>> fe_values_curvature;
  std::vector<Tensor<1, dim>>    projected_phase_indicator_gradient_values;
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
  // Used to evaluate the void fraction at particle locations (arbitrary
  // points within the cell); stored as a shared_ptr because it is only
  // instantiated when needed
  std::shared_ptr<FEPointEvaluation<1, dim>> point_evaluator_void_fraction;

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
  // Store the particle-location curls with the aliases above so the call to
  // FEValuesViews::Vector::get_function_curls() receives the exact curl type
  // expected by the installed deal.II version.
  std::vector<velocity_curl_type_2d>
    fluid_velocity_curls_at_particle_location_2d;
  std::vector<velocity_curl_type_3d>
                          fluid_velocity_curls_at_particle_location_3d;
  std::vector<Point<dim>> particle_reference_location;
  std::vector<double>     particle_weights;
  std::vector<double>     cell_void_fraction;
  std::vector<double>     Re_particle;
  unsigned int            max_number_of_particles_per_cell;
  unsigned int            number_of_particles;
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
   * Scratch component for the ALE method
   */
  bool gather_ale;

  /**
   * Scratch components for the mortar method
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
