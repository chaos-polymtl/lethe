// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
   * @brief Setup the degrees of freedom
   *
   */
  void
  setup_dofs() override;


  /**
   * @brief Establishes the constraints of the void fraction systems.
   *
   * @param boundary_conditions The boundary conditions of fluid dynamics. This is used to establish periodic boundary conditions for the void fraction.
   *
   */
  void
  setup_constraints(
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions);


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
   * @brief Initializes the void fraction at the beggining of a simulation
   *
   * @param time The current time. This is used for time-dependent functions.
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

  /// Index set for the locally owned degree of freedoms
  IndexSet locally_owned_dofs;

  /// Index set for the locally relevant degree of freedoms
  IndexSet locally_relevant_dofs;

  /// Constraints used for the boundary conditions of the void fraction.
  /// Currently, this is only used to establish periodic void fractions. This
  /// object has to be made public because the boundary conditions are set
  /// outside of the object for now.
  AffineConstraints<double> void_fraction_constraints;


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
   * @param is_post_mesh_adaptation Indicates if the equation is being solved during post_mesh_adaptation(),
   * for verbosity.
   */
  virtual void
  solve_linear_system_and_update_solution(
    const bool &is_post_mesh_adaptation = false) override;

  /**
   * @brief This function calculates and returns the periodic offset distance of the domain which is needed
   * for the periodic boundary conditions using the QCM or SPM for void fraction
   * with the GLS VANS/CFD-DEM solver. The distance is based on one of the
   * periodic boundaries and all particle location shifted by this distance is
   * according to this periodic boundary.
   *
   * @param boundary_id The id of one of the periodic boundaries
   *
   * @return The periodic offset distance
   */
  inline Tensor<1, dim>
  get_periodic_offset_distance(unsigned int boundary_id) const
  {
    Tensor<1, dim> offset;

    // Iterating over the active cells in the triangulation
    for (const auto &cell : (*this->triangulation).active_cell_iterators())
      {
        if (cell->is_locally_owned() || cell->is_ghost())
          {
            if (cell->at_boundary())
              {
                // Iterating over cell faces
                for (unsigned int face_id = 0; face_id < cell->n_faces();
                     ++face_id)
                  {
                    unsigned int face_boundary_id =
                      cell->face(face_id)->boundary_id();

                    // Check if face is on the boundary, if so, get
                    // the periodic offset distance for one pair of periodic
                    // faces only since periodic boundaries are aligned with the
                    // direction and only axis are currently allowed
                    if (face_boundary_id == boundary_id)
                      {
                        Point<dim> face_center = cell->face(face_id)->center();
                        auto periodic_cell = cell->periodic_neighbor(face_id);
                        unsigned int periodic_face_id =
                          cell->periodic_neighbor_face_no(face_id);
                        Point<dim> periodic_face_center =
                          periodic_cell->face(periodic_face_id)->center();

                        offset = periodic_face_center - face_center;

                        return offset;
                      }
                  }
              }
          }
      }

    // A zero tensor is returned in case no cells are found on the periodic
    // boundaries on this processor. This processor won't handle particle in
    // cells at periodic boundaries, so it won't affect any computation.
    return offset;
  }

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

  /// Locally owned solution of the void fraction
  GlobalVectorType void_fraction_locally_owned;

  /// Preconditioner used for the solution of the smoothed L2 projection of the
  /// void fraction
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /// System matrix used to assemble the smoothed L2 projection of the void
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

  /// Boolean to indicate if the mesh has periodic boundaries
  bool has_periodic_boundaries;

  /// Offset for the periodic boundary condition
  Tensor<1, dim> periodic_offset;

  /// Direction associated with the periodic boundary condition
  unsigned int periodic_direction;

  /// Vertices to periodic cell map, this is used in the QCM and the satellite
  /// point method to access the neighboring cells of a cell.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;
};



#endif
