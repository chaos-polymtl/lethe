/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */


#ifndef lethe_void_fraction_h
#define lethe_void_fraction_h

#include <core/vector.h>

#include <solvers/physics_subequations_solver.h>

#include <fem-dem/parameters_cfd_dem.h>

#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Calculates the area of intersection between a circular (2D) particle and a circle
 *
 * @param r_particle Radius of the particle
 *
 * @param r_circle Radius of the circle
 *
 * @param neighbor_distance Distance between the particle and the circle
 */
inline double
particle_circle_intersection_2d(double r_particle,
                                double r_circle,
                                double neighbor_distance)
{
  return pow(r_particle, 2) * Utilities::fixed_power<-1, double>(
                                cos((pow(neighbor_distance, 2) +
                                     pow(r_particle, 2) - pow(r_circle, 2)) /
                                    (2 * neighbor_distance * r_particle))) +
         Utilities::fixed_power<2, double>(r_circle) *
           Utilities::fixed_power<-1, double>(
             cos((pow(neighbor_distance, 2) - pow(r_particle, 2) +
                  pow(r_circle, 2)) /
                 (2 * neighbor_distance * r_circle))) -
         0.5 * sqrt((-neighbor_distance + r_particle + r_circle) *
                    (neighbor_distance + r_particle - r_circle) *
                    (neighbor_distance - r_particle + r_circle) *
                    (neighbor_distance + r_particle + r_circle));
}

/**
 * @brief Calculates the volume of intersection between a spherical (3D) particle and a sphere
 *
 * @param r_particle Radius of the particle
 *
 * @param r_sphere Radius of the sphere
 *
 * @param neighbor_distance Distance between the particle and the sphere
 */

inline double
particle_sphere_intersection_3d(double r_particle,
                                double r_sphere,
                                double neighbor_distance)
{
  return M_PI *
         Utilities::fixed_power<2, double>(r_sphere + r_particle -
                                           neighbor_distance) *
         (Utilities::fixed_power<2, double>(neighbor_distance) +
          (2 * neighbor_distance * r_particle) -
          (3 * Utilities::fixed_power<2, double>(r_particle)) +
          (2 * neighbor_distance * r_sphere) + (6 * r_sphere * r_particle) -
          (3 * Utilities::fixed_power<2, double>(r_sphere))) /
         (12 * neighbor_distance);
}

/**
 * @brief Void fraction calculator. This class manages the calculation of the
 * void fraction which is used as an auxiliary field for the solvers that solve
 * the Volume-Averaged Navier-Stokes equations. The present architecture
 * of the class support all of the different void fraction calculation methods
 * within a single class instead of building a class hierarchy. This is
 * because there are only a few void fraction calculation strategies.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class VoidFractionBase : public PhysicsLinearSubequationsSolver
{
public:
  VoidFractionBase(
    parallel::DistributedTriangulationBase<dim>             *triangulation,
    std::shared_ptr<Parameters::VoidFractionParameters<dim>> input_parameters,
    const Parameters::LinearSolver  &linear_solver_parameters,
    Particles::ParticleHandler<dim> *particle_handler,
    const unsigned int               fe_degree,
    const bool                       simplex,
    const ConditionalOStream        &pcout)
    : PhysicsLinearSubequationsSolver(pcout)
    , dof_handler(*triangulation)
    , triangulation(triangulation)
    , void_fraction_parameters(input_parameters)
    , linear_solver_parameters(linear_solver_parameters)
    // BB TODO verify if the particle handler can be remade const
    , particle_handler(particle_handler)
  {
    if (simplex)
      {
        fe         = std::make_shared<FE_SimplexP<dim>>(fe_degree);
        mapping    = std::make_shared<MappingFE<dim>>(*fe);
        quadrature = std::make_shared<QGaussSimplex<dim>>(fe->degree + 1);
      }
    else
      {
        // Usual case, for quad/hex meshes
        fe         = std::make_shared<FE_Q<dim>>(fe_degree);
        mapping    = std::make_shared<MappingQ<dim>>(fe->degree);
        quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
      }
  }

  /**
   * @brief Setup the degrees of freedom and establishes the constraints of the void fraction systems.
   *
   */
  void
  setup_dofs() override;


  /**
   * @brief Assembles the diagonal of the mass matrix to impose constraints on the system
   *
   *  @param diagonal_mass_matrix The matrix for which the diagonal entries will be filled with the diagonal of the mass matrix
   *
   *  @todo Establish if it is really necessary to keep this as a matrix and not as a vector since a diagonal is nothing more than a vector.
   */
  void
  assemble_mass_matrix_diagonal(
    TrilinosWrappers::SparseMatrix &diagonal_mass_matrix);

  /**
   * @brief Calculates the void fraction
   *
   * @param time current time for which the void fraction is to be calculated.
   *
   */
  void
  calculate_void_fraction(const double time);


  /**
   * @brief Assemble and solve the system.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  void
  solve(const bool & /*is_post_mesh_adaptation*/) override
  {}


  /**
   * @brief Percolates the time vector for the void fraction. This operation is called at the end of a time-step.
   *
   */
  void
  percolate_void_fraction()
  {
    for (unsigned int i = previous_void_fraction.size() - 1; i > 0; --i)
      {
        previous_void_fraction[i] = previous_void_fraction[i - 1];
      }
    previous_void_fraction[0] = void_fraction_locally_relevant;
  }


  /**
   * @brief Initializes the void fraction at the begginig of a simulation
   *
   */

  void
  initialize_void_fraction(const double time)
  {
    calculate_void_fraction(time);
    for (auto &previous_solution : this->previous_void_fraction)
      previous_solution = void_fraction_locally_relevant;
  }

  /// DoFHandler that manages the void fraction
  DoFHandler<dim> dof_handler;

  /// The solutions are made public instead of using getters
  /// Solution of the void fraction at previous time steps
  std::vector<GlobalVectorType> previous_void_fraction;

  /// Fully distributed (including locally relevant) solution
  GlobalVectorType void_fraction_locally_relevant;

  /// Finite element for the void fraction
  std::shared_ptr<FiniteElement<dim>> fe;

private:
  /**
   * @brief Calculates the void fraction using a function. This is a straightforward usage of VectorTools.
   *
   * @param time current time for which the void fraction is to be calculated.
   *
   */
  void
  calculate_void_fraction_function(const double time);

  /**
   * @brief Calculates the void fraction using the particle centered method.
   *
   */
  void
  calculate_void_fraction_particle_centered_method();

  /**
   * @brief Calculates the void fraction using the satellite point method.
   *
   */
  void
  calculate_void_fraction_satellite_point_method();

  /**
   * @brief Calculates the void fraction using the Quadrature-Centered Method (QCM).
   *
   */
  void
  calculate_void_fraction_quadrature_centered_method();

  /**
   * @brief Solve the linear system resulting from the assemblies
   *
   */
  virtual void
  solve_linear_system_and_update_solution(
    const bool &is_post_mesh_adaptation = false) override;

  /// Triangulation
  parallel::DistributedTriangulationBase<dim> *triangulation;

  /// Mapping for the void fraction
  std::shared_ptr<Mapping<dim>> mapping;

  /// Quadrature for the void fraction
  std::shared_ptr<Quadrature<dim>> quadrature;

  /// Parameters for the calculation of the void fraction
  std::shared_ptr<Parameters::VoidFractionParameters<dim>>
    void_fraction_parameters;

  /// Linear solvers for the calculation of the void fraction
  const Parameters::LinearSolver linear_solver_parameters;

  /// Particle handler used when the void fraction depends on particles
  Particles::ParticleHandler<dim> *particle_handler;



  /// Index set for the locally owned degree of freedoms
  IndexSet locally_owned_dofs;

  /// Index set for the locally relevant degree of freedoms
  IndexSet locally_relevant_dofs;



  /// Locally owned solution of the void fraction
  GlobalVectorType void_fraction_locally_owned;

  /// ??? Not fully sure of this yet
  TrilinosWrappers::SparseMatrix complete_system_matrix_void_fraction;
  GlobalVectorType               complete_system_rhs_void_fraction;

  /// Mass matrix used to constraint the value of the void fraction to be
  /// bounded
  TrilinosWrappers::SparseMatrix mass_matrix;

  /// Mass matrix diagonal used to constraint the value of the void fraction to
  /// be bounded
  GlobalVectorType diagonal_of_mass_matrix;

  /// ??? BB Not fully sure of this yet
  IndexSet active_set;

  /// Preconditioner used for the solution of the smoothed L2 projection of the
  /// void fraction
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /// Constraints used for the boundary conditions of the void fraction.
  /// Currently, this is only used to establish periodic void fractions.
  AffineConstraints<double> void_fraction_constraints;

  /// System matrix used to assembled the smoothed L2 projection of the void
  /// fraction
  TrilinosWrappers::SparseMatrix system_matrix_void_fraction;

  /// Right-hand side used to assemble the smoothed L2 projection of the void
  /// fraction
  GlobalVectorType system_rhs_void_fraction;

  /// Vertices to cell map, this is used in the QCM and the satellite point
  /// method to access the neighboring cells of a cell.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;


  /**
   * Member related to boundary conditions. At the moment, only a single
   * boundary condition is supported.
   */

  /// Boolean to indicate if the mesh has periodic boundaries
  bool has_periodic_boundaries;

  /// Offset for the periodic boundary condition
  Tensor<1, dim> periodic_offset;

  /// Direction associated of the periodic boundary condition
  unsigned int periodic_direction;

  /// Vertices to periodic cell map, this is used in the QCM and the satellite
  /// point method to access the neighboring cells of a cell.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;
};



#endif
