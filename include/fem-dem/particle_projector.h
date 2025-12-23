// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_void_fraction_h
#define lethe_void_fraction_h

#include <core/parameters_cfd_dem.h>
#include <core/vector.h>

#include <solvers/navier_stokes_scratch_data.h>
#include <solvers/physics_subequations_solver.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Calculate the area of intersection between a circular (2D) particle and a circle.
 *
 * @param[in] r_particle Radius of the particle
 *
 * @param[in] r_circle Radius of the circle
 *
 * @param[in] neighbor_distance Distance between the particle and the circle
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
 * @brief Calculate the volume of intersection between a spherical (3D) particle and a sphere.
 *
 * @param[in] r_particle Radius of the particle
 * @param[in] r_sphere Radius of the sphere
 * @param[in] neighbor_distance Distance between the particle and the sphere
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
 * @brief ParticleFieldQCM calculator.
 * This class stores the required information for the calculation of a
 * field or a particle-fluid interaction onto a mesh. This class does not solve
 * any equation, but compartmentalize the necessary information for the
 * projection of the particle field or interaction onto a mesh. Multiple
 * instances of this function can be called to establish new projections.In the
 * case of a particle-fluid interaction, there is no particle field to project
 * and, consequently, the property_start_index parameter is given an invalid
 * default value (-1).
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam n_components The number of components in the field. This number should be either 1 (a scalar) or dim (a Tensor<1,dim>)
 * @tparam property_start_index An integer that indicates at which particle property the field starts.
 * If the propery_start_index is negative, then there is no property start and
 * this project is not used to project a particle field.
 */
template <int dim, int n_components, int property_start_index = -1>
class ParticleFieldQCM
{
public:
  /**
   * @brief Constructor
   *
   * @param triangulation The triangulation on which the particles reside
   * @param fe_degree The finite element degree used to interpolate the field
   * @param simplex A flag to indicate if the simulations are being done with simplex elements.
   * @param neumann_boundaries A flag to indicate what to do with cells that do not contain particles.
   * If the neumann_boundaries is set to true, then  cells without particles do
   * not have a mass matrix assembled. This corresponds to setting a neumann
   * boundary at the cells which do not contain particles.
   * @param conservative_projection A flag to indicate if the projection
   * of the field must be conservative. When a projection is conservative,
   * the integral of the field over the triangulation will be strictly
   * equal to the field summed over all particles. This is useful when
   * projecting a quantity such as the force acting on the particles. When
   * conservative is put to false, the projection is meant to be a smooth
   * representation of the particle information and does not conserve it.
   */
  ParticleFieldQCM(parallel::DistributedTriangulationBase<dim> *triangulation,
                   const unsigned int                           fe_degree,
                   const bool                                   simplex,
                   const bool neumann_boundaries,
                   const bool conservative_projection)
    : dof_handler(*triangulation)
    , neumann_boundaries(neumann_boundaries)
    , conservative_projection(conservative_projection)
  {
    if (simplex)
      {
        const FE_SimplexP<dim> fe_temp(fe_degree);
        fe = std::make_shared<FESystem<dim>>(fe_temp, n_components);
      }
    else
      {
        const FE_Q<dim> temp_fe(fe_degree);
        fe = std::make_shared<FESystem<dim>>(temp_fe, n_components);
      }
  }

  /// DoFHandler that manages the void fraction
  DoFHandler<dim> dof_handler;

  /// deal.II vector for the particle velocity
  LinearAlgebra::distributed::Vector<double> particle_field_solution;

  /// Finite element for the particle field
  std::shared_ptr<FESystem<dim>> fe;

  /// Index set for the locally owned degree of freedoms
  IndexSet locally_owned_dofs;

  /// Index set for the locally relevant degree of freedoms
  IndexSet locally_relevant_dofs;

  /// System matrix used to assemble the smoothed L2 projection of the void
  /// fraction
  TrilinosWrappers::SparseMatrix system_matrix;

  /// Right-hand side used to assemble the smoothed L2 projection of the void
  /// fraction
  LinearAlgebra::distributed::Vector<double> system_rhs;

  /// Constraints used for the boundary conditions of the particle field
  AffineConstraints<double> particle_field_constraints;

  /// Boolean that indicates if Neumann boundary conditions are used during the
  /// projection In this case, instead of assuming a value of zero if there are
  /// no particles, the mass matrix is cancelled in region without particles
  /// which results in a no flux zone.
  const bool neumann_boundaries;

  const bool conservative_projection;

  /**
   * @brief Setup the degrees of freedom. This function allocates the necessary memory.
   */
  void
  setup_dofs();
};

/**
 * @brief Particle information projection.
 * This class manages the calculation of the projection of the particle
 * fields onto a triangulation. Its main use is the calculation of the
 * void fraction which is used as an auxiliary field for the solvers that solve
 * the Volume-Averaged Navier-Stokes equations. The present architecture
 * of the class support all of the different void fraction calculation methods
 * within a single class instead of building a class hierarchy. This is
 * because there are only a few void fraction calculation strategies.
 *
 * Using QCM, this class is also capable of projecting arbitrary fields
 * belonging to the particle back onto the mesh by instantiating
 * @ParticleFieldQCM class.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class ParticleProjector : public PhysicsLinearSubequationsSolver
{
public:
  /**
   * @brief Constructor for the void fraction calculator.
   *
   * @param triangulation The triangulation on which the simulation is being done.
   * @param input_parameters The parameters for the void fraction calculation.
   * @param linear_solver_parameters The parameters for the linear solver. This is used to solve the L2 projection of the void fraction
   * @param particle_handler The particle handler used to manage the particles in the simulation. This object must be provided even if the particles are not used to establish the void fraction.
   * @param fe_degree The finite element degree used to interpolate the void fraction.
   * @param simplex A flag to indicate if the simulations are being done with simplex elements.
   * @param pcout The ConditionalOStream used to print the information.
   */
  ParticleProjector(
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
    , particle_handler(particle_handler)
    , particle_have_been_projected(false)
    , particle_velocity(triangulation, fe_degree, simplex, true, false)
    , fluid_force_on_particles_two_way_coupling(triangulation,
                                                fe_degree,
                                                simplex,
                                                false,
                                                true)
    , fluid_drag_on_particles(triangulation, fe_degree, simplex, false, true)
    , momentum_transfer_coefficient(triangulation,
                                    fe_degree,
                                    simplex,
                                    false,
                                    true)
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
        fe      = std::make_shared<FE_Q<dim>>(fe_degree);
        mapping = std::make_shared<MappingQ<dim>>(fe->degree);
        if (this->void_fraction_parameters->quadrature_rule ==
            Parameters::VoidFractionQuadratureRule::gauss)
          {
            if (this->void_fraction_parameters->n_quadrature_points == 0)
              quadrature = std::make_shared<QGauss<dim>>(fe->degree + 1);
            else
              quadrature = std::make_shared<QGauss<dim>>(
                this->void_fraction_parameters->n_quadrature_points);
          }
        else if (this->void_fraction_parameters->quadrature_rule ==
                 Parameters::VoidFractionQuadratureRule::gauss_lobatto)
          {
            if (this->void_fraction_parameters->n_quadrature_points == 0)
              quadrature = std::make_shared<QGaussLobatto<dim>>(fe->degree + 2);
            else if (this->void_fraction_parameters->n_quadrature_points >= 3)
              quadrature = std::make_shared<QGaussLobatto<dim>>(
                this->void_fraction_parameters->n_quadrature_points);
            else
              throw(std::runtime_error(
                "For void fraction using Gauss-Lobatto ('gauss-lobatto') quadrature rule, the minimum number of quadrature points is 3"));
          }
        else
          {
            AssertThrow(
              false,
              ExcMessage(
                "An invalid quadrature rule was provided to the ParticleProject class."));
          }
      }

    AssertThrow(
      !void_fraction_parameters->project_particle_velocity ||
        void_fraction_parameters->l2_smoothing_length > 0,
      ExcMessage(
        "The projection of the particle velocity field requires that the l2 smoothing length be > 0. Otherwise, the Poisson equation used in the Neumann boundary is invalid."));

    AssertThrow(
      void_fraction_parameters->project_particle_velocity == false ||
        (void_fraction_parameters->project_particle_velocity &&
         void_fraction_parameters->mode == Parameters::VoidFractionMode::qcm),
      ExcMessage(
        "The projection of the particle velocity currently requires that the QCM method be used for the calculation of the void fraction."));
  }

  /**
   * @brief Setup the degrees of freedom.
   *
   */
  void
  setup_dofs() override;

  /**
   * @brief Establish the constraints of the void fraction systems.
   *
   * @param[in] boundary_conditions The boundary conditions of fluid dynamics.
   * This is used to establish periodic boundary conditions for the void
   * fraction.
   *
   */
  void
  setup_constraints(
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions);

  /**
   * @brief Calculate the void fraction.
   *
   * @param[in] time Current time for which the void fraction is to be
   * calculated.
   *
   */
  void
  calculate_void_fraction(const double time);

  /**
   * @brief Calculate all the necessary information for the particle-fluid coupling.
   * This means calculating the projection of the particle-fluid interaction
   * forces as well as the solid velocity onto the CFD mesh.
   *
   */
  template <typename VectorType>
  void
  calculate_particle_fluid_forces_projection(
    const Parameters::CFDDEM      &cfd_dem_parameters,
    DoFHandler<dim>               &fluid_dof_handler,
    const VectorType              &present_velocity_pressure_solution,
    const std::vector<VectorType> &previous_velocity_pressure_solution,
    const Tensor<1, 3>            &gravity,
    NavierStokesScratchData<dim>   scratch_data);

  /**
   * @brief Assemble and solve the system.
   */
  void
  solve() override
  {}

  /**
   * @brief Percolate the time vector for the void fraction. This operation is called at the end of a time-step. It percolates both the deal.II and the Trilinos vectors.
   *
   * TODO - Refactor the ParticleProjector class to use deal.II vectors for
   * everything instead of a blend of Trilinos and deal.II vector.
   *
   */
  void
  percolate_void_fraction()
  {
    for (unsigned int i = previous_void_fraction.size() - 1; i > 0; --i)
      {
        previous_void_fraction[i] = previous_void_fraction[i - 1];
        void_fraction_previous_solution[i] =
          void_fraction_previous_solution[i - 1];
      }
    previous_void_fraction[0]          = void_fraction_locally_relevant;
    void_fraction_previous_solution[0] = void_fraction_solution;
  }

  /**
   * @brief Initialize the void fraction at the beginning of a simulation.
   *
   * @param[in] time The current time. This is used for time-dependent
   * functions.
   */
  void
  initialize_void_fraction(const double time)
  {
    calculate_void_fraction(time);
    for (auto &previous_solution : this->previous_void_fraction)
      previous_solution = void_fraction_locally_relevant;

    for (auto &previous_solution : this->void_fraction_previous_solution)
      previous_solution = void_fraction_solution;
  }

  /**
   * @brief Calculates the projection of a particle field onto the mesh.
   *
   * @tparam property_start_index An integer for the particle_property that will be projected
   *
   * @tparam n_components The number of components of the field. This should be either 1 (scalar) or dim (a Tensor<1,dim>).
   *
   * @param field_qcm The container for the field that will be projected.
   *
   */

  template <int property_start_index, int n_components>
  void
  calculate_field_projection(
    ParticleFieldQCM<dim, property_start_index, n_components> &field_qcm);


  /**
   * @brief Zeros out and update the ghost values of the auxiliary fields.
   * This is necessary when the void fraction calculation mode does not rely on
   * particles.
   */
  void
  zero_out_and_ghost_auxiliary_fields()
  {
    fluid_force_on_particles_two_way_coupling.particle_field_solution = 0.;
    fluid_force_on_particles_two_way_coupling.particle_field_solution
      .update_ghost_values();
    momentum_transfer_coefficient.particle_field_solution = 0.;
    momentum_transfer_coefficient.particle_field_solution.update_ghost_values();
    particle_velocity.particle_field_solution = 0.;
    particle_velocity.particle_field_solution.update_ghost_values();
    fluid_drag_on_particles.particle_field_solution = 0.;
    fluid_drag_on_particles.particle_field_solution.update_ghost_values();
  }

  /// DoFHandler that manages the void fraction
  DoFHandler<dim> dof_handler;

  /// Mapping for the void fraction
  std::shared_ptr<Mapping<dim>> mapping;

  /// Quadrature for the void fraction
  std::shared_ptr<Quadrature<dim>> quadrature;

  /// The solutions are made public instead of using getters
  /// Solution of the void fraction at previous time steps
  std::vector<GlobalVectorType> previous_void_fraction;

  /// Fully distributed (including locally relevant) solution
  GlobalVectorType void_fraction_locally_relevant;

  /// Locally owned solution of the void fraction
  GlobalVectorType void_fraction_locally_owned;

  /// deal.II vector for the void fraction
  LinearAlgebra::distributed::Vector<double> void_fraction_solution;

  /// deal.II vector the previous void fraction
  std::vector<LinearAlgebra::distributed::Vector<double>>
    void_fraction_previous_solution;

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
   * @brief Calculate the characteristic radius of a sphere with a given volume (area in 2D):
   * R = (2*dim*V/pi)^(1/dim) / 2
   *
   * @param volume value for which the equivalent radius of a sphere is calculated.
   *
   * @return Radius of a sphere with equivalent volume
   */
  static double
  radius_sphere_equivalent_volume(const double volume)
  {
    return 0.5 * pow(2.0 * static_cast<double>(dim) * volume / M_PI,
                     1.0 / static_cast<double>(dim));
  }

  /**
   * @brief Calculate the radius of the QCM averaging sphere
   *
   * @param cell_measure The measure of the cell in wich QCM is calculated.
   *
   * @return The QCM radius used in the calculations.
   */
  inline double
  calculate_qcm_radius_from_cell_measure(const double cell_measure)
  {
    if (void_fraction_parameters->qcm_sphere_equal_cell_volume == true)
      {
        // Get the radius by the volume of sphere which is
        // equal to the volume of cell
        return radius_sphere_equivalent_volume(cell_measure);
      }
    else
      {
        // The radius is obtained from the volume of sphere based
        // on R_s = h_omega
        return std::cbrt(cell_measure);
      }
  }

  /**
   * @brief Calculate the intersection measure (volume in 3D and area in 2D) between a particle and an hypersphere.
   * This is used when calculating the intersection volume between the averaging
   * sphere and a particle in QCM. The calculation inherently assumes that the
   * radius of the averaging sphere is larger than that of the particle.
   *
   * @param r_particle The radius of the particle
   * @param r_sphere The radius of the averaging sphere
   * @param distance_between_spheres The distance between the centers of the spheres
   *
   * @return The intersection volume (in 3D) or area (in 2d)
   */

  static double
  calculate_intersection_measure(double r_particle,
                                 double r_sphere,
                                 double distance_between_spheres)
  {
    // Particle completely in the reference sphere
    if (distance_between_spheres <= (r_sphere - r_particle))
      {
        if constexpr (dim == 2)
          return M_PI * Utilities::fixed_power<2>(r_particle);
        if constexpr (dim == 3)
          return 4. / 3. * M_PI * Utilities::fixed_power<3>(r_particle);
      }

    // Particle partially in the reference sphere
    if ((distance_between_spheres > (r_sphere - r_particle)) &&
        (distance_between_spheres < (r_sphere + r_particle)))
      {
        if constexpr (dim == 2)
          return particle_circle_intersection_2d(r_particle,
                                                 r_sphere,
                                                 distance_between_spheres);

        else if constexpr (dim == 3)
          return particle_sphere_intersection_3d(r_particle,
                                                 r_sphere,
                                                 distance_between_spheres);
      }

    // Particle completely outside the reference
    // sphere. The intersection volume is zero.
    return 0;
  }

  /**
   * @brief Calculate the void fraction using a function. This is a straightforward usage of VectorTools.
   *
   * @param[in] time Current time for which the void fraction is to be
   * calculated.
   *
   */
  void
  calculate_void_fraction_function(const double time);

  /**
   * @brief Calculate the void fraction using the particle centered method.
   *
   */
  void
  calculate_void_fraction_particle_centered_method();

  /**
   * @brief Calculate the void fraction using the satellite point method.
   *
   */
  void
  calculate_void_fraction_satellite_point_method();

  /**
   * @brief Calculate the void fraction using the Quadrature-Centered Method (QCM).
   *
   */
  void
  calculate_void_fraction_quadrature_centered_method();

  /**
   * @brief Solve the linear system resulting from the assemblies.
   *
   */
  virtual void
  solve_linear_system_and_update_solution() override;

  /**
   * @brief Calculate and return the periodic offset distance vector of the domain which is needed
   * for the periodic boundary conditions using the QCM or SPM for void fraction
   * with the GLS VANS/CFD-DEM solver. The distance is based on one of the
   * periodic boundaries. This periodic boundary is then used as the reference
   * for all particle locations shifted by this offset vector.
   *
   * @param[in] boundary_id The ID of one of the periodic boundaries
   *
   * @return The periodic offset vector.
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

  /// Parameters for the calculation of the void fraction
  /// Right now this is used for the VANS matrix-free solver
  /// to directly calculate the void fraction from the function itself.
  /// The function will be removed back to private in a future PR.
  std::shared_ptr<Parameters::VoidFractionParameters<dim>>
    void_fraction_parameters;

  /// Linear solvers for the calculation of the void fraction
  const Parameters::LinearSolver linear_solver_parameters;

  /// Particle handler used when the void fraction depends on particles
  Particles::ParticleHandler<dim> *particle_handler;


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

  /// Smoothing length factor for the void fraction calculation
  const double l2_smoothing_factor =
    Utilities::fixed_power<2>(void_fraction_parameters->l2_smoothing_length);

  /// Boolean indicator used to check if the particles have at least been
  /// projected once. This is mainly used for assertions and sanity checking
  /// purposes
  bool particle_have_been_projected;

public:
  /// Projector used to save the particle velocity field onto the mesh
  /// This field is not volumetric, meaning that the integral over the
  /// entire domain of this field will not yield the sum of the particle
  /// velocity over all the particles.
  ParticleFieldQCM<dim, 3, DEM::CFDDEMProperties::PropertiesIndex::v_x>
    particle_velocity;

  /// Projector used to save the forces used in the two-way coupling scheme,
  /// without drag. This field is volumetric. The integral of this field over
  /// the entire domain will result in the sum of the force (without drag) over
  /// all of the particles.
  ParticleFieldQCM<
    dim,
    3,
    DEM::CFDDEMProperties::PropertiesIndex::fem_force_two_way_coupling_x>
    fluid_force_on_particles_two_way_coupling;

  /// Projector used to store the drag force. The drag is stored as a separate
  /// field since we may want to make implicit calculations of it. This field
  /// is volumetric. The integral of this field over the entire domain will
  /// result in the sum of the drag force over all the particles.
  ParticleFieldQCM<dim, 3, DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x>
    fluid_drag_on_particles;

  /// Projector used to save the momentum transfer coefficient (beta). This
  /// field is volumetric. The integral of this field over the entire domain
  /// will result in the sum of the momentum transfer coefficient over all
  /// particles.
  ParticleFieldQCM<
    dim,
    1,
    DEM::CFDDEMProperties::PropertiesIndex::momentum_transfer_coefficient>
    momentum_transfer_coefficient;
};
#endif
