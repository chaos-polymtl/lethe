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
///   [] Gather output hook
///   [] Finish simulation
///   [x] Percolate_time_vector
///   - modify_solution
///   [] update_boundary_conditions
///   [] get_dof_handler
///   [] postprocess
///   [] pre_mesh_adaptation
///   [] post_mesh_adaptation
///   - write_checkpoint
///   - read_checkpoint
///   - gather_tables()
///   - compute_kelly
///   - set_initial_conditions
///   [] setup_preconditioner
/// [] Physics field
/// [] Multiphysics interface components
/// - Physical Properties
///   - Electrical conductivity
///   - Permetivity
///   - Permeability
/// [] FEM section for DPG
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
    std::shared_ptr<SimulationControl> p_simulation_control)
    : AuxiliaryPhysics<dim, GlobalVectorType>()
    , multiphysics(multiphysics_interface)
    , computing_timer(p_triangulation->get_mpi_communicator(),
                      this->pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , dof_handler_trial_interior(
        std::make_shared<DoFHandler<dim>>(*triangulation))
    , dof_handler_trial_skeleton(
        std::make_shared<DoFHandler<dim>>(*triangulation))
    , dof_handler_test(std::make_shared<DoFHandler<dim>>(*triangulation))
    , extractor_E_real(0)
    , extractor_E_imag(dim)
    , extractor_H_real(2 * dim)
    , extractor_H_imag(3 * dim)
  {
    if (simulation_parameters.mesh.simplex)
      {
        // for simplex meshes
        AssertThrow(
          false,
          "TimeHarmonicMaxwell solver not yet implemented for simplex meshes.");
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe_trial_interior = std::make_shared<FESystem<dim>>(
          FE_DGQ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order) ^
            dim,
          FE_DGQ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order) ^
            dim,
          FE_DGQ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order) ^
            dim,
          FE_DGQ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order) ^
            dim);
        fe_trial_skeleton = std::make_shared<FESystem<dim>>(
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order));
        fe_test = std::make_shared<FESystem<dim>>(
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order),
          FE_NedelecSZ<dim>(
            simulation_parameters.fem_parameters.electromagnetic_order));
        mapping = std::make_shared<MappingQ<dim>>(fe_trial_interior->degree);
        cell_quadrature = std::make_shared<QGauss<dim>>(fe_test->degree + 1);
        face_quadrature =
          std::make_shared<QGauss<dim - 1>>(fe_test->degree + 1);
      }

    // Initialize solution shared_ptr
    present_solution = std::make_shared<GlobalVectorType>();

    // Allocate solution transfer
    solution_transfer =
      std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
        *dof_handler_trial_interior);
  }

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
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation()
  {}

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  virtual void
  percolate_time_vectors()
  {}

  /**
   * @brief Carry out modifications on the auxiliary physic solution.
   * To be defined for some physics only (eg. free surface, see vof.h).
   */
  virtual void
  modify_solution() {};

  /**
   * @brief Update non zero constraints if the boundary is time-dependent
   */
  virtual void
  update_boundary_conditions() {};

  /**
   * @brief Getter method to access the private attribute dof_handler for the
   * physic currently solved. NB : The dof_handler that is returned is the one
   * for the interior trial space only has it is where the solution lives. The
   * other two DoFHandlers of the class are used for computation only.
   */
  virtual const DoFHandler<dim> &
  get_dof_handler()
  {
    return *dof_handler_trial_interior;
  }

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess(bool first_iteration) {};


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  virtual void
  pre_mesh_adaptation() {};

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  virtual void
  post_mesh_adaptation() {};

  /**
   * @brief Prepares auxiliary physics to write checkpoint
   */
  virtual void
  write_checkpoint() {};

  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  virtual void
  read_checkpoint() {};

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for a given TimeHarmonicMaxwell solver.
   *
   * @return Structure containing the TableHandler objects and their corresponding file names.
   */
  virtual std::vector<OutputStructTableHandler>
  gather_tables()
  {
    return std::vector<OutputStructTableHandler>();
  };

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
                dealii::Vector<float> &estimated_error_per_cell) {};

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() {};

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionally support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() {};

  /**
   * @brief Set up preconditioner. Not used for the auxiliary physics but
   * needed for the compilation of the non-linear solver.
   */
  void
  setup_preconditioner() override {};

private:
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
   * @brief A vector containing all the values of the solution.
   */
  std::shared_ptr<GlobalVectorType> present_solution;

  /**
   * @brief The right hand side vector.
   */
  GlobalVectorType system_rhs;

  /**
   * @brief Store the nonzero constraints that arise from several sources such
   * as boundary conditions and hanging nodes in the mesh. See the deal.II
   * documentation on constraints on degrees of freedom for more information.
   */
  AffineConstraints<double> constraints;

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
