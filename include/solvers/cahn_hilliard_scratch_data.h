// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cahn_hilliard_scratch_data_h
#define lethe_cahn_hilliard_scratch_data_h


#include <solvers/multiphysics_interface.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/fe/fe_system.h>

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
class CahnHilliardScratchData : public PhysicsScratchDataBase
{
public:
  /**
   * @brief Constructor. The constructor creates the fe_values that will be used
   * to fill the member variables of the scratch. It also allocated the
   * necessary memory for all member variables. However, it does not do any
   * evaluation, since this needs to be done at the cell level.
   *
   * @param properties_manager The physical properties manager (see physical_properties_manager.h)
   * @param fe_cahn_hilliard The FESystem used to solve the Cahn-Hilliard equations
   * @param quadrature The quadrature to use for the assembly
   * @param mapping The mapping of the domain in which the Cahn-Hilliard equations are solved
   * @param fe_fd The FESystem used to solve the Fluid Dynamics equations
   *
   */
  CahnHilliardScratchData(const PhysicalPropertiesManager &properties_manager,
                          const FESystem<dim>             &fe_cahn_hilliard,
                          const Quadrature<dim>           &quadrature,
                          const Mapping<dim>              &mapping,
                          const FiniteElement<dim>        &fe_fd,
                          const Quadrature<dim - 1>       &face_quadrature);

  /**
   * @brief Copy Constructor. Same as the main constructor.
   * This constructor only uses the other scratch to build the FeValues, it
   * does not copy the content of the other scratch into itself since, by
   * definition of the WorkStream mechanism it is assumed that the content of
   * the scratch will be reset on a cell basis.
   *
   * @param sd
   */
  CahnHilliardScratchData(const CahnHilliardScratchData<dim> &sd);


  /** @brief Allocates the necessary memory for all members of the scratch
   */
  void
  allocate() override;

  /** @brief Reinitializes the content of the scratch.
   *
   * Using the FeValues and the content of the solutions and previous solutions,
   * fills all of the class member of the scratch.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the fe which is used to fill the
   * FeValues.
   *
   * @param[in] current_solution The present value of the solution for [Phi,eta]
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   * @param[in] source_function The function describing the source term in
   * Cahn-Hilliard equations.
   *
   */

  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         Function<dim>                 *source_function);

  /**
   * @brief Reinitialize the velocity, calculated by the Fluid Dynamics.
   *
   * @tparam VectorType The Vector type used for the solvers.
   *
   * @param[in] cell The cell for which the velocity is reinitialized. This cell
   * must be compatible with the Fluid Dynamics FE.
   *
   * @param[in] current_solution The present value of the solution for
   * \f$[u,p]\f$.
   *
   * @param[in] ale ALE parameters.
   */
  template <typename VectorType>
  void
  reinit_velocity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const VectorType           &current_solution,
                  const Parameters::ALE<dim> &ale);

  /**
   * @brief Calculates the physical properties that may be required by the
   * Cahn-Hilliard problem. Namely the surface tension and the Cahn-Hilliard
   * mobility.
   */
  void
  calculate_physical_properties();


  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  dealii::types::material_id           material_id;
  std::vector<double>                  density;
  std::vector<double>                  kinematic_viscosity;
  std::vector<double>                  surface_tension;
  std::vector<double>                  mobility_cahn_hilliard;
  std::vector<double>                  mobility_cahn_hilliard_gradient;


  FEValuesExtractors::Scalar phase_order;
  FEValuesExtractors::Scalar chemical_potential;

  // FEValues for the Cahn-Hilliard problem
  FEValues<dim> fe_values_cahn_hilliard;
  unsigned int  n_dofs;
  unsigned int  n_q_points;
  double        cell_size;

  // Quadrature
  Table<1, double>        JxW;
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
  Table<2, double>         phi_phase;
  Table<2, Tensor<2, dim>> hess_phi_phase;
  Table<2, double>         laplacian_phi_phase;
  Table<2, Tensor<1, dim>> grad_phi_phase;
  Table<2, double>         phi_potential;
  Table<2, Tensor<2, dim>> hess_phi_potential;
  Table<2, double>         laplacian_phi_potential;
  Table<2, Tensor<1, dim>> grad_phi_potential;

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
  FEFaceValues<dim> fe_face_values_cahn_hilliard;
  Table<2, double>  face_JxW;

  unsigned int n_faces;
  unsigned int n_faces_q_points;

  // Is boundary cell indicator
  bool                      is_boundary_cell;
  std::vector<bool>         is_boundary_face;
  std::vector<unsigned int> boundary_face_id;

  // First is face number, second quadrature point, third DOF
  Table<3, Tensor<1, dim>> grad_phi_face_phase;
  // First vector is face number, second quadrature point
  std::vector<std::vector<Tensor<1, dim>>> face_phase_grad_values;
  // The normal vector is necessary for the free angle boundary condition
  std::vector<std::vector<Tensor<1, dim>>> face_normal;
};

#endif
