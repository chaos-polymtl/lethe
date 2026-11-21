// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_scratch_data_h
#define lethe_tracer_scratch_data_h

#include <core/ale.h>
#include <core/lethe_grid_tools.h>

#include <solvers/physical_properties_manager.h>
#include <solvers/physics_scratch_data.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_system.h>

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
class TracerScratchData : public PhysicsScratchDataBase
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
   * @param uses_dg If true, the scratch is used for the discontinuous Galerkin
   * formulation. The hessians, laplacians and previous gradients, which are
   * only required by the stabilization of the continuous formulation, are then
   * not computed.
   *
   */
  TracerScratchData(const PhysicalPropertiesManager &properties_manager,
                    const FiniteElement<dim>        &fe_tracer,
                    const Quadrature<dim>           &quadrature,
                    const Quadrature<dim - 1>       &face_quadrature,
                    const Mapping<dim>              &mapping,
                    const FiniteElement<dim>        &fe_fd,
                    const bool                       uses_dg);

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
  TracerScratchData(const TracerScratchData<dim> &sd);


  /** @brief Allocates the memory for the scratch
   *
   * This function allocates the necessary memory for all members of the scratch
   *
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
   * @param[in] current_solution The present value of the solution.
   *
   * @param[in] previous_solutions The solutions at the previous time steps.
   *
   * @param[in] source_function The function describing the tracer source term.
   */
  template <typename VectorType>
  void
  reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
         const VectorType                                     &current_solution,
         const std::vector<VectorType> &previous_solutions,
         const Function<dim>           *source_function);

  /**
   * @brief Reinitializes the signed distance vector of the scratch.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * @param[in] immersed_solid_shape Shared pointer to the shape describing the
   * particles (if there are any).
   */
  void
  reinit_signed_distance(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const std::shared_ptr<Shape<dim>>                    &immersed_solid_shape);

  /** @brief Reinitializes the signed distance vector of the scratch for a given face.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * @param[in] face_no The face index associated with the cell
   * @param[in] immersed_solid_shape Shared pointer to the shape describing the
   * particles (if there are any).
   */
  void
  reinit_signed_distance_at_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                   &face_no,
    const std::shared_ptr<Shape<dim>>                    &immersed_solid_shape);


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
    const typename DoFHandler<dim>::active_cell_iterator  &cell,
    const VectorType                                      &velocity_solution,
    const Parameters::ALE<dim>                            &ale,
    const std::shared_ptr<Functions::ParsedFunction<dim>> &drift_velocity);


  /** @brief Reinitializes the content of the scratch for the internal faces. This is only used for the DG assemblers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face_no The face index associated with the cell
   *
   * @param[in] sub_face_no The subface index associated with the face
   *
   * @param[in] neigh_cell The neighboring cell
   *
   * @param[in] neigh_face_no The face index associated with the neighboring
   * cell
   *
   * @param[in] neigh_sub_face_no The subface index associated with the
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
    const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
    const unsigned int                                   &neigh_face_no,
    const unsigned int                                   &neigh_sub_face_no,
    const VectorType                                     &current_solution);


  /** @brief Reinitializes the content of the scratch for the boundary faces. This is only used for the DG assemblers.
   *
   * @param[in] cell The cell over which the assembly is being carried.
   *
   * @param[in] face The face index associated with the cell
   *
   * @param[in] current_solution The present value of the solution.
   * there are any).
   */
  template <typename VectorType>
  void
  reinit_boundary_face(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const unsigned int                                   &face_no,
    const unsigned int                                   &boundary_index,
    const VectorType                                     &current_solution);

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
    const typename DoFHandler<dim>::active_cell_iterator  &velocity_cell,
    const unsigned int                                    &face_no,
    const VectorType                                      &velocity_solution,
    const Parameters::ALE<dim>                            &ale,
    const std::shared_ptr<Functions::ParsedFunction<dim>> &drift_velocity);

  /**
   * @brief Renitialize the content of the scratch data for mortar
   *
   * @param[in] cell The cell over which the assembly is being carried.
   * This cell must be compatible with the FE which is used to fill the
   * FeValues.
   *
   * @param[in] mortar_parameters Parameters for the mortar method
   *
   * @param[in] time Current simulation time
   */
  void
  reinit_mortar(const typename DoFHandler<dim>::active_cell_iterator &cell,
                const Parameters::Mortar<dim> &mortar_parameters,
                const double                  &time);

  /** @brief Calculates the physical properties at a face
   */
  void
  calculate_face_physical_properties();


  /** @brief Calculates the physical properties. This function calculates the physical properties
   * that may be required by the tracer. Namely the diffusivity.
   *
   */
  void
  calculate_physical_properties();

  // True if the scratch is used for the discontinuous Galerkin formulation
  const bool uses_dg;

  // Physical properties
  PhysicalPropertiesManager            properties_manager;
  std::map<field, std::vector<double>> fields;
  std::vector<double>                  tracer_diffusivity;
  std::vector<double>                  tracer_diffusivity_face;
  std::vector<double>                  tracer_diffusivity_0;
  std::vector<double>                  tracer_diffusivity_1;
  std::vector<double>                  tracer_reaction_prefactor;

  // Gradient of the tracer reaction prefactor with respect to the concentration
  // This is calculated by deriving the prefactor by the concentration
  // (dk/dC)
  std::vector<double> grad_tracer_reaction_prefactor;

  // FEValues for the Tracer problem
  FEValues<dim>          fe_values_tracer;
  FEInterfaceValues<dim> fe_interface_values_tracer;

  unsigned int n_dofs;
  unsigned int n_interface_dofs;
  unsigned int n_q_points;
  double       cell_size;

  // Quadrature
  Table<1, double>        JxW;
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

  // Shape function quantities at the internal faces, indexed by quadrature
  // point and interface dof
  Table<2, double> face_phi_here;
  Table<2, double> face_phi_there;
  Table<2, double> face_jump_phi;
  Table<2, double> face_average_grad_phi_dot_n;

  // SIPG (interior faces) or Nitsche (boundary) penalization factor
  // The penalty factor is calculated using a measure of the element size. It is
  // multiplied by the diffusion coefficient to apply the SIPG method
  double penalty_factor;

  // Bool that defines if the selected face is a dirichlet/outlet boundary
  unsigned int boundary_index;

  // Immersed solid shape (signed distance function)
  std::vector<double> sdf_values;

  // Source term
  std::vector<double> source;

  // Shape functions
  Table<2, double>         phi;
  Table<2, Tensor<2, dim>> hess_phi;
  Table<2, double>         laplacian_phi;
  Table<2, Tensor<1, dim>> grad_phi;

  /**
   * Scratch components for the mortar method
   */
  bool                        gather_mortar;
  std::vector<Tensor<1, dim>> rotor_linear_velocity_values;

  /**
   * Scratch component for the Navier-Stokes component
   */
  FEValuesExtractors::Vector velocities;
  // This FEValues must mandatorily be instantiated for the velocity
  FEValues<dim>               fe_values_fd;
  FEFaceValues<dim>           fe_face_values_fd;
  std::vector<Tensor<1, dim>> velocity_values;
  std::vector<double>         velocity_divergences;
  std::vector<Tensor<1, dim>> face_velocity_values;
};

#endif
