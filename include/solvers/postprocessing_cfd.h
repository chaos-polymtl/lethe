/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Audrey Collard-Daigneault, Bruno Blais, Polytechnique Montreal, 2020
 -
 */


#ifndef lethe_postprocessing_cfd_h


// Base
#  include <deal.II/base/quadrature_lib.h>

// Lac
#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#  include <deal.II/lac/vector.h>

// Dofs
#  include <deal.II/dofs/dof_handler.h>

// Fe
#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/mapping_fe.h>

// Lethe includes
#  include <core/boundary_conditions.h>
#  include <core/parameters.h>

#  include <solvers/physical_properties_manager.h>

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
calculate_pressure_drop(const DoFHandler<dim> &       dof_handler,
                        std::shared_ptr<Mapping<dim>> mapping,
                        const VectorType &            evaluation_point,
                        const Quadrature<dim> &       cell_quadrature_formula,
                        const Quadrature<dim - 1> &   face_quadrature_formula,
                        const unsigned int            inlet_boundary_id,
                        const unsigned int            outlet_boundary_id);

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
              const VectorType &     evaluation_point,
              const double           time_step,
              const Quadrature<dim> &quadrature_formula,
              const Mapping<dim> &   mapping);

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
                    const VectorType &     evaluation_point,
                    const Quadrature<dim> &quadrature_formula,
                    const Mapping<dim> &   mapping);

/**
 * @brief Calculate the average kinetic energy in the simulation domain
 * @return average kinetic energy in the domain
 * Post-processing function
 * This function calculates the average kinetic energy in the simulation domain
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
calculate_kinetic_energy(const DoFHandler<dim> &dof_handler,
                         const VectorType &     evaluation_point,
                         const Quadrature<dim> &quadrature_formula,
                         const Mapping<dim> &   mapping);

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
calculate_apparent_viscosity(const DoFHandler<dim> &    dof_handler,
                             const VectorType &         evaluation_point,
                             const Quadrature<dim> &    quadrature_formula,
                             const Mapping<dim> &       mapping,
                             PhysicalPropertiesManager &properties_manager);

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
std::vector<std::vector<Tensor<1, dim>>>
calculate_forces(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  PhysicalPropertiesManager &                          properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1> &                          face_quadrature_formula,
  const Mapping<dim> &                                 mapping);


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
std::vector<Tensor<1, 3>>
calculate_torques(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  PhysicalPropertiesManager &                          properties_manager,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const Quadrature<dim - 1> &                          face_quadrature_formula,
  const Mapping<dim> &                                 mapping);


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
                   const VectorType &     evaluation_point,
                   const Function<dim> *  exact_solution,
                   const Quadrature<dim> &quadrature_formula,
                   const Mapping<dim> &   mapping);


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
calculate_flow_rate(const DoFHandler<dim> &    dof_handler,
                    const VectorType &         present_solution,
                    const unsigned int &       boundary_id,
                    const Quadrature<dim - 1> &face_quadrature_formula,
                    const Mapping<dim> &       mapping);



#  define lethe_postprocessing_cfd_h

#endif
