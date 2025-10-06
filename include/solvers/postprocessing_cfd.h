// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessing_cfd_h
#define lethe_postprocessing_cfd_h

#include <core/boundary_conditions.h>

#include <solvers/physical_properties_manager.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

/**
 * @brief Calculate the pressure drop and total pressure drop between two boundaries.
 * The pressure drop thus calculated has units of Length^2/Time^2.
 *
 * @return Pressure drop of the flow between two boundaries of the domain
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution for which the pressure drop is calculated. The velocity field is assumed to be the "dim" field
 *
 * @param cell_quadrature_formula The cell quadrature formula for the calculation
 *
 * @param face_quadrature_formula The face quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 *
 * @param inlet_boundary_id The id of the inlet boundary
 *
 * @param outlet_boundary_id The id of the outlet boundary
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_pressure_drop(const DoFHandler<dim>     &dof_handler,
                        const Mapping<dim>        &mapping,
                        const VectorType          &evaluation_point,
                        const Quadrature<dim>     &cell_quadrature_formula,
                        const Quadrature<dim - 1> &face_quadrature_formula,
                        const unsigned int         inlet_boundary_id,
                        const unsigned int         outlet_boundary_id);

/**
 * @brief Calculate the CFL condition on the simulation domain
 * @return CFL maximal value in the domain
 * Post-processing function
 * This function calculates the maximal CFL value in the domain
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution for which the CFL is calculated. The velocity field is assumed to be the first field.
 *
 * @param quadrature_formula The quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_CFL(const DoFHandler<dim> &dof_handler,
              const VectorType      &evaluation_point,
              const double           time_step,
              const Quadrature<dim> &quadrature_formula,
              const Mapping<dim>    &mapping);

/**
 * @brief Calculate the average enstrophy in the simulation domain
 * @return average kinetic energy in the domain
 * Post-processing function
 * This function calculates the average enstrophy in the simulation domain
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution at which the force is calculated
 *
 * @param quadrature_formula The quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_enstrophy(const DoFHandler<dim> &dof_handler,
                    const VectorType      &evaluation_point,
                    const Quadrature<dim> &quadrature_formula,
                    const Mapping<dim>    &mapping);

/**
 * @brief Calculate the average kinetic energy in the simulation domain
 * @return average kinetic energy in the domain
 * Post-processing function
 * This function calculates the average kinetic energy in the simulation domain
 *
 * @param[in] dof_handler The dof_handler used for the calculation
 *
 * @param[in] evaluation_point The solution for the calculation
 *
 * @param[in] quadrature_formula The quadrature formula for the calculation
 *
 * @param[in] mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_kinetic_energy(const DoFHandler<dim> &dof_handler,
                         const VectorType      &evaluation_point,
                         const Quadrature<dim> &quadrature_formula,
                         const Mapping<dim>    &mapping);

/**
 * @brief Calculate the average power done by pressure in the simulation
 * domain. The average power is defined as ∫u.∇pdΩ/∫1dΩ
 * @return Average power done by pressure in the simulation domain
 * Post-processing function
 *
 * @param[in] dof_handler The dof_handler used for the calculation
 *
 * @param[in] evaluation_point The solution for the calculation
 *
 * @param[in] quadrature_formula The quadrature formula for the calculation
 *
 * @param[in] mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_pressure_power(const DoFHandler<dim> &dof_handler,
                         const VectorType      &evaluation_point,
                         const Quadrature<dim> &quadrature_formula,
                         const Mapping<dim>    &mapping);

/**
 * @brief Calculate the viscous dissipation of kinetic energy which is defined as
 *  ∫∇u.τdΩ/∫1dΩ
 * @return Viscous dissipation of kinetic energy
 * Post-processing function
 *
 * @param[in] dof_handler The dof_handler used for the calculation
 *
 * @param[in] evaluation_point The solution for the calculation
 *
 * @param[in] quadrature_formula The quadrature formula for the calculation
 *
 * @param[in] mapping The mapping of the simulation
 *
 * @param[in] properties_manager Manager for the physical properties used to
 * calculate the kinematic viscosity
 */
template <int dim, typename VectorType>
double
calculate_viscous_dissipation(
  const DoFHandler<dim>           &dof_handler,
  const VectorType                &evaluation_point,
  const Quadrature<dim>           &quadrature_formula,
  const Mapping<dim>              &mapping,
  const PhysicalPropertiesManager &properties_manager);

/**
 * @brief Calculates the apparent viscosity of the fluid for non Newtonian flows.
 * @return the apparent viscosity
 * Post-processing frunction
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution at which the force is calculated
 *
 * @param quadrature_formula The quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 *
 * @param physical_properties The parameters containing the required physical properties
 */
template <int dim, typename VectorType>
double
calculate_apparent_viscosity(
  const DoFHandler<dim>           &dof_handler,
  const VectorType                &evaluation_point,
  const Quadrature<dim>           &quadrature_formula,
  const Mapping<dim>              &mapping,
  const PhysicalPropertiesManager &properties_manager);

/**
 * @brief Calculates the force due to the fluid motion on every boundary conditions
 * @return std::vector of forces on each boundary condition
 * Post-processing function
 * This function calculates the force acting on each of the boundary conditions
 * within the domain. It generates a vector which size is the number of boundary
 * conditions
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution at which the force is calculated
 *
 * @param physical_properties The parameters containing the required physical properties
 *
 * @param boundary_conditions The boundary conditions object
 *
 * @param face_quadrature_formula The face quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
std::vector<std::map<types::boundary_id, Tensor<1, dim>>>
calculate_forces(
  const DoFHandler<dim>                               &dof_handler,
  const VectorType                                    &evaluation_point,
  const PhysicalPropertiesManager                     &properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1>                           &face_quadrature_formula,
  const Mapping<dim>                                  &mapping);


/**
 * @brief Calculates the torques due to the fluid motion on every boundary conditions
 * @return std::vector of torques on each boundary condition
 * Post-processing function
 * This function calculates the torqueacting on each of the boundary conditions
 * within the domain. It generates a vector which size is the number of boundary
 * conditions.
 *
 * @param dof_handler The dof_handler used for the calculation.
 *
 * @param evaluation_point The solution at which the torque is calculated.
 *
 * @param physical_properties The parameters containing the required physical properties.
 *
 * @param boundary_conditions The boundary conditions object.
 *
 * @param face_quadrature_formula The face quadrature formula for the calculation.
 *
 * @param mapping The mapping of the simulation.
 */
template <int dim, typename VectorType>
std::map<types::boundary_id, Tensor<1, 3>>
calculate_torques(
  const DoFHandler<dim>                               &dof_handler,
  const VectorType                                    &evaluation_point,
  const PhysicalPropertiesManager                     &properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1>                           &face_quadrature_formula,
  const Mapping<dim>                                  &mapping);


/**
 * @brief Calculates the L2 norm of the error on velocity and pressure
 * @return std::pair<double,double> containing the L2 norm of the error for velocity and pressure
 * Post-processing function
 * This function calculates the L2 norm of the error on the velocity and
 * pressure. Since the pressure in GLS is defined up to a given constant, the
 * mean solution of the pressure is removed from both the analytical and the
 * numerical solution to ensure that the convergence can be monitored correctly.
 *
 * @param dof_handler The dof_handler used for the calculation.
 *
 * @param evaluation_point The solution at which the torque is calculated.
 *
 * @param exact_solution The exact solution, a function of dim+1 component for velocity + pressure
 *
 * @param quadrature_formula The quadrature formula for the calculation.
 *
 * @param mapping The mapping of the simulation.
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_L2_error(const DoFHandler<dim> &dof_handler,
                   const VectorType      &evaluation_point,
                   const Function<dim>   *exact_solution,
                   const Quadrature<dim> &quadrature_formula,
                   const Mapping<dim>    &mapping);

/**
 * @brief calculate_flow_rate. This function calculates the volumetric flow
 * rate at the selected boundary. It actually calculates the flow rate through
 * the summation of the value at each cell surface with the normal vector,
 * the velocity value and the area.
 *
 * @param dof_handler. The argument used to get finite elements
 *
 * @param present_solution. The vector which contains all the values to
 *                          calculate the flow rate
 *
 * @param boundary_id. The inlet boundary
 *
 * @param face_quadrature_formula The face quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
std::pair<double, double>
calculate_flow_rate(const DoFHandler<dim>     &dof_handler,
                    const VectorType          &present_solution,
                    const unsigned int        &boundary_id,
                    const Quadrature<dim - 1> &face_quadrature_formula,
                    const Mapping<dim>        &mapping);

/**
 * @brief calculate_average_velocity. This function calculates average velocity
 * from the volumetric flow rate at the selected boundary.
 *
 * @param dof_handler. The argument used to get finite elements
 *
 * @param present_solution. The vector which contains all the values to
 *                          calculate the flow rate
 *
 * @param boundary_id. The inlet boundary
 *
 * @param face_quadrature_formula The face quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_average_velocity(const DoFHandler<dim>     &dof_handler,
                           const VectorType          &present_solution,
                           const unsigned int        &boundary_id,
                           const Quadrature<dim - 1> &face_quadrature_formula,
                           const Mapping<dim>        &mapping);

/**
 * @brief calculate_average_velocity. This function calculates the average velocity of
 * the domain for the solver GLS VANS.
 *
 * Ū = ∫ε(v⋅γ)dΩ/∫εdΩ
 * Where Ū is the average velocity, ε is the void fraction, v is the velocity
 * of the cell, γ is the vector of the cell (flow direction) and dΩ is the
 * volume of the cell.
 *
 * @param dof_handler. The argument used to get velocity at quadrature points
 *
 * @param void_fraction_dof_handler. The argument used to get void fraction at quadrature points
 *
 * @param present_solution. The vector which contains the velocity values
 *
 * @param present_void_fraction_solution. The vector which contains the void fraction values
 *
 * @param flow_direction. The flow direction
 *
 * @param quadrature_formula The quadrature formula for the calculation
 *
 * @param mapping The mapping of the simulation
 */
template <int dim, typename VectorType>
double
calculate_average_velocity(const DoFHandler<dim> &dof_handler,
                           const DoFHandler<dim> &void_fraction_dof_handler,
                           const VectorType      &present_solution,
                           const VectorType   &present_void_fraction_solution,
                           const unsigned int &flow_direction,
                           const Quadrature<dim> &quadrature_formula,
                           const Mapping<dim>    &mapping);



#define lethe_postprocessing_cfd_h

#endif
