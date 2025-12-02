// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_h
#define lethe_time_harmonic_maxwell_h

#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/data_out.h>

#include <memory.h>



/**
 * @class TimeHamonicMaxwell Auxiliary physics solver for the time harmonic
 * form of the Maxwell equations.
 * TODO Add more details here.
 *
 **/

/// TODO in this PR:
/// [] Remplir chacune des fonctions
///   [x] Constructeur
///   [x] Destructeur
///   [] Setup_dofs
///   [x] Gather output hook
///   [x] Calculate_L2_error
///   [x] Finish simulation
///   [x] Percolate_time_vector
///   - modify_solution
///   [] update_boundary_conditions
///   [] get_dof_handler
///   [x] postprocess
///   [] pre_mesh_adaptation
///   [] post_mesh_adaptation
///   - write_checkpoint
///   - read_checkpoint
///   - gather_tables()
///   - compute_kelly
///   - set_initial_conditions
///   [] setup_preconditioner
/// [x] Physics field
/// [x] Multiphysics interface components
/// [] Physical Properties
///   [] Electrical conductivity
///   [] Permetivity
///   [] Permeability
/// [X] FEM section for DPG
/// - Mesh adaptation

using VectorType = GlobalVectorType;
template <int dim>
class TimeHarmonicMaxwell : public AuxiliaryPhysics<dim, VectorType>
{
public:
  /**
   * @brief Constructor for the TimeHarmonicMaxwell object
   *
   * @param multiphysics_interface Map of the auxiliary physics that will be
   * solved on top of a computational fluid dynamic simulation.
   *
   * @param p_simulation_parameters Contain the simulation parameter file
   * information.
   *
   * @param p_triangulation Contain the mesh information. In a
   * parallel::DistributedTriangulationBase<dim> not every detail may be known
   * on each processor. The mesh is distributed between the processors.
   *
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
   * @brief Carry out modifications on the auxiliary physic solution.
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
   * @brief Compute the Kelly error estimator used to refine mesh on a auxiliary physic parameter.
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
   * @brief Compute the DPG error estimator based on the energy norm residual used to refine mesh of the TimeHarmonicMaxwell physic.
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
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionally support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() override
  {
    AssertThrow(
      false,
      ExcMessage(
        " The TimeHarmonicMaxwell solver does not support initial conditions as it is a frequency domain solver. Consequently, this method always solves a steady-state problem."));
  };

  /**
   * @brief Set up preconditioner. Not used for the auxiliary physics but
   * needed for the compilation of the non-linear solver.
   */
  void
  setup_preconditioner() override;

  /**
   * @brief Define the constraints for the TimeHarmonicMaxwell solver.
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
   * physic currently solved. NB : The dof_handler that is returned is the one
   * for the interior trial space only has it is where the solution lives. The
   * other two DoFHandlers of the class are used for computation only.
   */
  virtual const DoFHandler<dim> &
  get_dof_handler() override
  {
    return *dof_handler_trial_interior;
  }

  /**
   * @brief Getter method to access the private attribute evaluation_point for
   * the physic currently solved.
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
   * newton_update for the physic currently solved.
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
   * present_solution for the physic currently solved. NB : present_solution is
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
   * present_solution_skeleton for the physic currently solved. NB :
   * present_solution_skeleton is not passed to the multiphysics interface at
   * the end of the setup_dofs method. In the case of TimeHarmonicMaxwell, it
   * corresponds to the solution on the skeleton for the trial space.
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
   * system_rhs for the physic currently solved.
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
   * nonzero_constraints for the physic currently solved.
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
  output_newton_update_norms(const unsigned int display_precision) override
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
   *  @brief Assemble the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

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
   * cells of the triangulation of the trial space for the interior dofs.
   */
  std::shared_ptr<DoFHandler<dim>> dof_handler_trial_interior;

  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation of the trial space for the skeleton dofs.
   */
  std::shared_ptr<DoFHandler<dim>> dof_handler_trial_skeleton;

  /**
   * @brief Given a triangulation and a description of a finite element, this
   * class enumerates degrees of freedom on all vertices, edges, faces, and
   * cells of the triangulation of the test space dofs. Note that in DPG, the
   * test space is not split into interior and skeleton dofs in opposition to
   * the trial space.
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
   * @brief Approximate an integral by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim>> cell_quadrature;
  /**
   * @brief Approximate an integral by evaluating the integrand at specific
   * points and summing the point values with specific weights.
   */
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;


  /////////////               Solution storage                /////////////

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
