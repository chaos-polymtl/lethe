// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_h
#define lethe_time_harmonic_maxwell_h

#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <memory.h>

#include <complex>
#include <numbers>


/**
 * @class TimeHarmonicMaxwell Auxiliary physics solver for the time harmonic
 * form of the Maxwell equations.
 * This class implements an auxiliary physics solver for the ultraweak
 * formulation of the time-harmonic Maxwell equations using the Discontinuous
 * Petrov-Galerkin (DPG) method. The solver is designed to work in
 * three-dimensional spaces and requires the 3 sets of Finite Element spaces
 * that define the DPG method: the interior trial space that is discontinuous
 * between elements, the skeleton trace space that lives on the faces of the
 *mesh elements and connects the interior solutions of adjacent elements, and
 *the test space that is enriched compared to the trial space and is also
 * discontinuous between elements.
 **/

/// TODO in this PR:
/// [] TimeHarmonicMaxwell class
///   [x] Constructeur
///   [x] Destructeur
///   [x] Gather output hook
///   [x] Calculate_L2_error
///   [x] Finish simulation
///   [x] Percolate_time_vector
///   [x] modify_solution
///   [x] update_boundary_conditions
///   [x] postprocess
///   [x] pre_mesh_adaptation
///   [x] post_mesh_adaptation
///   - write_checkpoint
///   - read_checkpoint
///   - gather_tables()
///   - compute_kelly
///   - compute_energy_norm
///   [x] Setup_dofs
///   [x] set_initial_conditions
///   [x] setup_preconditioner
///   [x] define_constraints
///   [x] solve_linear_system
///   [x] get_dof_handler
///   [x] get_evaluation_point
///   [x] get_local_evaluation_point
///   [x] get_newton_update
///   [x] get_present_solution
///   [x] get_present_solution_skeleton
///   [x] get_system_rhs
///   [x] get_nonzero_constraints
///   [x] output_newton_update_norms
///   [x] get_residual_rescale_metric
///   [x] assemble_system_matrix
///   [x] assemble_system_rhs
///   - assemble_local_system_matrix
///   - assemble_local_system_rhs
///   - setup_assemblers
///   - copy_local_matrix_to_global_matrix
///   - copy_local_rhs_to_global_rhs
/// [x] Physics field
/// [x] Multiphysics interface components
/// - Physical Properties
///   - Electrical conductivity
///   - Permittivity
///   - Permeability
/// [X] FEM section for DPG
/// - Mesh adaptation

using VectorType = GlobalVectorType;

DeclException1(
  TimeHarmonicMaxwellBoundaryConditionMissing,
  types::boundary_id,
  << "The boundary id: " << arg1
  << " is defined in the triangulation, but not as a boundary condition for the TimeHarmonicMaxwell physics. Lethe does not assign a default boundary condition to boundary ids. Every boundary id defined within the triangulation must have a corresponding boundary condition defined in the input file.");

DeclException1(
  TimeHarmonicMaxwellDimensionNotSupported,
  int,
  << "The time-harmonic Maxwell solver does not support dimension: " << arg1
  << ". Currently, only 3D problems are supported as the 2D version of curls and cross products have completely different definitions than their 3D counterparts.");

template <int dim>
class TimeHarmonicMaxwell : public AuxiliaryPhysics<dim, VectorType>
{
public:
  /**
   * @brief Constructor for the TimeHarmonicMaxwell object
   *
   * @param multiphysics_interface Map of the auxiliary physics that will be
   * solved on top of a computational fluid dynamic simulation.
   * @param p_simulation_parameters Contains the simulation parameter file
   * information.
   * @param p_triangulation Contains the mesh information. In a
   * parallel::DistributedTriangulationBase<dim> not every detail may be known
   * on each processor. The mesh is distributed between the processors.
   * @param p_simulation_control Object responsible for the control of
   * steady-state and transient simulations. Contains all the information
   * related to time stepping and the stopping criteria.
   *
   */
  TimeHarmonicMaxwell(
    MultiphysicsInterface<dim>      *multiphysics_interface,
    const SimulationParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control);

  /**
   * @brief TimeHarmonicMaxwell - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  virtual ~TimeHarmonicMaxwell()
  {}

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  virtual std::vector<OutputStruct<dim, VectorType>>
  gather_output_hook() override;

  /**
   * @brief Calculates the L2 error of the solution. It returns a vector of doubles containing the L2 errors of the real and imaginary parts of the electric and magnetic fields in the order [E_real, E_imag, H_real, H_imag].
   */
  std::vector<double>
  calculate_L2_error();

  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation() override;

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  virtual void
  percolate_time_vectors() override;

  /**
   * @brief Carry out modifications on the auxiliary physics solution.
   * To be defined for some physics only (eg. free surface, see vof.h).
   */
  virtual void
  modify_solution() override;

  /**
   * @brief Update non zero constraints if the boundary is time-dependent
   */
  virtual void
  update_boundary_conditions() override;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess(bool first_iteration) override;


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  virtual void
  pre_mesh_adaptation() override;

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  virtual void
  post_mesh_adaptation() override;

  /**
   * @brief Prepares auxiliary physics to write checkpoint
   */
  virtual void
  write_checkpoint() override;

  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  virtual void
  read_checkpoint() override;

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for a given TimeHarmonicMaxwell solver.
   *
   * @return Structure containing the TableHandler objects and their corresponding file names.
   */
  virtual std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Compute the Kelly error estimator used to refine mesh on an auxiliary physics parameter.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  virtual void
  compute_kelly(const std::pair<const Variable,
                                Parameters::MultipleAdaptationParameters> &ivar,
                dealii::Vector<float> &estimated_error_per_cell) override;

  /**
   * @brief Compute the DPG error estimator based on the energy norm residual used to refine mesh of the TimeHarmonicMaxwell physics.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  virtual void
  compute_energy_norm(
    const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                          &ivar,
    dealii::Vector<float> &estimated_error_per_cell);

  /**
   * @brief Sets up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Sets up the initial conditions associated with the physics. Generally, most physics
   * only support imposing nodal values, but some physics additionally support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() override;

  /**
   * @brief Sets up preconditioner of the system matrix.
   */
  void
  setup_preconditioner() override;

  /**
   * @brief Define the constraints for the TimeHarmonicMaxwell solver. This is equivalent to the define_non_zero_constraints method in other physics solvers, but here there is no need to separate zero and non-zero constraints.
   */
  void
  define_constraints();

  /**
   * @brief Call for the solution of the linear system of equation using CG solver.
   */
  void
  solve_linear_system() override;

  /**
   * @brief Getter method to access the private attribute dof_handler for the
   * physics currently solved. NB : The dof_handler that is returned by this
   * function is the one for the interior trial space of the DPG system. It is
   * where the interior solution lives and therefore is the one we are
   * interested in. The other two DoFHandlers of the class are used for
   * computation to obtain the interior of the cell solution.
   */
  virtual const DoFHandler<dim> &
  get_dof_handler() override
  {
    return *dof_handler_trial_interior;
  }

  /**
   * @brief Getter method to access the private attribute evaluation_point for
   * the currently solved physics.
   *
   * @return The vector at which the evaluation is performed.
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    AssertThrow(
      false,
      ExcMessage(
        "The TimeHarmonicMaxwell solver is linear, this method should not be called."));
  }

  /**
   * @brief Getter method to access the private attribute
   * local_evaluation_point for the physic currently solved.
   *
   * @return The local evaluation point. Ghosts cells are not considered in
   * this evaluation.
   */
  GlobalVectorType &
  get_local_evaluation_point() override
  {
    AssertThrow(
      false,
      ExcMessage(
        "The TimeHarmonicMaxwell solver is linear, this method should not be called."));
  }

  /**
   * @brief Getter method to access the private attribute
   * newton_update for the currently solved physics.
   *
   * @return The direction used to perform the newton iteration.
   */
  GlobalVectorType &
  get_newton_update() override
  {
    AssertThrow(
      false,
      ExcMessage(
        "The TimeHarmonicMaxwell solver is linear, this method should not be called."));
  }

  /**
   * @brief Getter method to access the private attribute
   * present_solution for the currently solved physics. NB : present_solution is
   * now passed to the multiphysics interface at the end of the setup_dofs
   * method. In the case of TimeHarmonicMaxwell, it corresponds to the solution
   * in the interior trial space.
   *
   * @return A vector containing all the values of the solution.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    return *present_solution;
  }


  /**
   * @brief Getter method to access the private attribute
   * present_solution_skeleton for the currently solved physics. NB :
   * present_solution_skeleton is only defined for the currently solved physics
   * and is not passed to the multiphysics interface at the end of the
   * setup_dofs method.
   *
   * @return A vector containing all the values of the solution.
   */
  GlobalVectorType &
  get_present_solution_skeleton()
  {
    return *present_solution_skeleton;
  }


  /**
   * @brief Getter method to access the private attribute
   * system_rhs for the currently solved physics.
   *
   * @return Right hand side vector.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    return system_rhs;
  }

  /**
   * @brief Getter method to access the private attribute
   * nonzero_constraints for the currently solved physics.
   *
   * @return The nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }

  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  void
  output_newton_update_norms(
    [[maybe_unused]] const unsigned int display_precision) override
  {
    AssertThrow(
      false,
      ExcMessage(
        "The TimeHarmonicMaxwell solver is linear, this method should not be called."));
  }

  /**
   * @brief Return the metric for residual rescaling. By default, should return 1.
   * If the rescale_residual_by_volume is set to true, the method
   * returns the global volume of the triangulation.
   *
   * @return Rescale metric.
   */
  double
  get_residual_rescale_metric() const override
  {
    return simulation_parameters.linear_solver.at(PhysicsID::VOF)
               .rescale_residual_by_volume ?
             std::sqrt(
               GridTools::volume(*this->triangulation, *this->mapping)) :
             1.;
  }

private:
  ///////                Physics solver core functions               ////////

  /**
   * @brief Verify consistency of the input parameters for boundary
   * conditions to ensure that for every boundary id within the
   * triangulation, a boundary condition has been specified in the input file.
   */
  void
  verify_consistency_of_boundary_conditions()
  {
    // Sanity check all of the boundary conditions of the triangulation to
    // ensure that they have a type.
    std::vector<types::boundary_id> boundary_ids_in_triangulation =
      this->triangulation->get_boundary_ids();
    for (auto const &boundary_id_in_tria : boundary_ids_in_triangulation)
      {
        AssertThrow(
          simulation_parameters
              .boundary_conditions_time_harmonic_electromagnetics.type.find(
                boundary_id_in_tria) !=
            simulation_parameters
              .boundary_conditions_time_harmonic_electromagnetics.type.end(),
          TimeHarmonicMaxwellBoundaryConditionMissing(boundary_id_in_tria));
      }
  }

  /**
   *  @brief Assemble the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Reconstruct the interior solution from the skeleton solution after solving
   * the linear system on the skeleton only.
   */
  void
  reconstruct_interior_solution();


  /**
   * This helper function projects a 3D tensor onto the tangential plane
   * defined by the given normal vector. Mathematically, it computes:
   * \f[
   * \mathbf{t} = \mathbf{n} \times (\mathbf{v} \times \mathbf{n})
   * \f]
   * which removes the normal component of the vector \f$\mathbf{v}\f$,
   * keeping only the tangential part.
   *
   * @tparam dim Spatial dimension.
   * @param tensor Input vector (field value) to be projected.
   * @param normal Unit normal vector defining the face orientation.
   * @return The tangential component of the input vector.
   *
   * @note This operation is used to obtain traces in H^{-1/2}(curl) spaces,
   * where only tangential components are continuous across interfaces. The
   * inline keyword is enforced to make sure that the tensor operations are
   * efficient as they are called frequently during assembly.
   */
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, std::complex<double>>
  map_H12(const Tensor<1, dim, std::complex<double>> &tensor,
          const Tensor<1, dim>                       &normal)
  {
    if (dim != 3)
      {
        AssertThrow(
          false,
          ExcMessage(
            "The map_H12 function is only implemented for 3D problems."));
      }
    auto result = cross_product_3d(normal, cross_product_3d(tensor, normal));

    return result;
  }

  /////// Auxiliary physics parameters for TimeHarmonicMaxwell solver ////////

  MultiphysicsInterface<dim> *multiphysics;

  /**
   * @brief Store information related to the computing time such as CPU times or
   * wall time.
   */
  TimerOutput computing_timer;

  /**
   * @brief Contain the simulation parameter file information.
   */
  const SimulationParameters<dim> &simulation_parameters;


  ////////////// Core elements of the TimeHarmonicMaxwell solver //////////////

  /**
   * @brief Collection of cells that cover the domain on which one wants to
   * solve a partial differential equation.
   */
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  /**
   * @brief Responsible for the control of steady-state and transient
   * simulations. Contains all the information related to time stepping and the
   * stopping criteria. See simulation_control abstract class for more
   * information.
   */
  std::shared_ptr<SimulationControl> simulation_control;

  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation of the trial space for the interior dofs. NB:
   * Since it is based on DG elements, there are no dofs on vertices, edges or
   * faces.
   */
  std::shared_ptr<DoFHandler<dim>> dof_handler_trial_interior;

  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation of the trial space for the skeleton dofs. NB:
   * This DoFHandler is based on global FE spaces with frozen interior dofs, so
   * the dofs effectively live only on faces, edges and vertices.
   */
  std::shared_ptr<DoFHandler<dim>> dof_handler_trial_skeleton;

  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation of the test space dofs. NB: in DPG, the
   * test space is not split into interior and skeleton dofs in opposition to
   * the trial space and therefore manages vertices, edges, faces and cells
   * dofs.
   */
  std::shared_ptr<DoFHandler<dim>> dof_handler_test;

  /**
   * @brief The base class for finite elements. This class will manage everything related to the shape functions of the trial interior space.
   */
  std::shared_ptr<FiniteElement<dim>> fe_trial_interior;

  /**
   * @brief The base class for finite elements. This class will manage everything related to the shape functions of the trial skeleton space.
   */
  std::shared_ptr<FiniteElement<dim>> fe_trial_skeleton;

  /**
   * @brief The base class for finite elements. This class will manage everything related to the shape functions of the test space.
   */
  std::shared_ptr<FiniteElement<dim>> fe_test;

  /**
   * @brief Preconditioner for the linear solver. At the present moment, this is an identity preconditioner.
   */
  std::shared_ptr<TrilinosWrappers::PreconditionIdentity> preconditioner;

  /**
   * @brief Store some convergence data, such as residuals of the cg-method,
   * or some evaluated <i>L<sup>2</sup></i>-errors of discrete solutions.
   * Evaluate convergence rates or orders.
   */
  ConvergenceTable error_table;

  /////////////               Mapping and Quadrature             /////////////

  /**
   * @brief Transformation which maps point in the reference cell to
   * points in the actual grid cell.
   */
  std::shared_ptr<Mapping<dim>> mapping;
  /**
   * @brief Approximate an integral over a cell by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim>> cell_quadrature;
  /**
   * @brief Approximate an integral over a cell's face by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;


  /////////////               Solution storage                /////////////

  /**
   * @brief IndexSet of the owned degrees of freedom for the interior trial space.
   */
  IndexSet locally_owned_dofs_trial_interior;

  /**
   * @brief IndexSet of the relevant degrees of freedom for the interior trial space.
   */
  IndexSet locally_relevant_dofs_trial_interior;

  /**
   * @brief IndexSet of the owned degrees of freedom for the skeleton trial space.
   */
  IndexSet locally_owned_dofs_trial_skeleton;

  /**
   * @brief IndexSet of the relevant degrees of freedom for the skeleton trial space.
   */
  IndexSet locally_relevant_dofs_trial_skeleton;

  /**
   * @brief IndexSet of the owned degrees of freedom for the test space.
   */
  IndexSet locally_owned_dofs_test;

  /**
   * @brief IndexSet of the relevant degrees of freedom for the test space.
   */
  IndexSet locally_relevant_dofs_test;

  /**
   * @brief The system matrix.
   */
  TrilinosWrappers::SparseMatrix system_matrix;

  /**
   * @brief A vector containing all the values of the solution in the interior trial space.
   */
  std::shared_ptr<GlobalVectorType> present_solution;

  /**
   * @brief A vector containing all the values of the solution on the skeleton for the trial space.
   */
  std::shared_ptr<GlobalVectorType> present_solution_skeleton;

  /*
   * @brief A vector containing all the values of the DPG built-in a-posteriori error indicator.
   */
  std::shared_ptr<GlobalVectorType> present_DPG_error_indicator;

  /**
   * @brief The right hand side vector.
   */
  GlobalVectorType system_rhs;

  /**
   * @brief Store the nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   * The zero constraints are not used in this solver as it solve a linear
   * system.
   */
  AffineConstraints<double> nonzero_constraints;

  /**
   * @brief SolutionTransfer<dim, GlobalVectorType>> is
   * used to implement the transfer of a discrete FE function
   * (e.g. a solution vector) from one mesh to another. This Deal.ii class is
   * used for mesh_refinement and simulation restarts.
   */
  std::shared_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;


  ///////         TimeHarmonicMaxwell specific parameters         ////////

  /**
   * @brief Extractor for the real part of the electric field vector.
   */
  const FEValuesExtractors::Vector extractor_E_real;
  /**
   * @brief Extractor for the imaginary part of the electric field vector.
   */
  const FEValuesExtractors::Vector extractor_E_imag;
  /**
   * @brief Extractor for the real part of the magnetic field vector.
   */
  const FEValuesExtractors::Vector extractor_H_real;
  /**
   * @brief Extractor for the imaginary part of the magnetic field vector.
   */
  const FEValuesExtractors::Vector extractor_H_imag;
};



#endif
