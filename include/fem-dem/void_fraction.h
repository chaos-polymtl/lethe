// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_void_fraction_h
#define lethe_void_fraction_h

#include <core/lethe_grid_tools.h>
#include <core/vector.h>

#include <solvers/physics_subequations_solver.h>

#include <fem-dem/parameters_cfd_dem.h>

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
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

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
 *
 * @param[in] r_sphere Radius of the sphere
 *
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
 * @brief Particle velocity calculator. This class stores the required information for the calculation of the projected particle velocity.
 * This class does not solve any equation, but compartimentalize the necessary
 * information for the projection of the particle velocity onto the mesh.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim, int component_start, int n_components>
class ParticleFieldQCM
{
public:
  ParticleFieldQCM(parallel::DistributedTriangulationBase<dim> *triangulation,
                   const unsigned int                           fe_degree,
                   const bool                                   simplex)
    : dof_handler(*triangulation)
  {
    if (simplex)
      {
        if constexpr (n_components == 1)
          fe = FE_SimplexP<dim>(fe_degree);
        else
          {
            const FE_SimplexP<dim> fe_temp(fe_degree);
            fe = std::make_shared<FESystem<dim>>(fe_temp, n_components);
          }
      }
    else
      {
        if constexpr (n_components == 1)
          fe = FE_Q<dim>(fe_degree);
        else
          {
            const FE_Q<dim> temp_fe(fe_degree);
            fe = std::make_shared<FESystem<dim>>(temp_fe, dim);
          }
      }
  }

  /// DoFHandler that manages the void fraction
  DoFHandler<dim> dof_handler;

  /// Fully distributed (including locally relevant) solution
  GlobalVectorType particle_field_locally_relevant;

  /// deal.II vector for the particle velocity
  LinearAlgebra::distributed::Vector<double> particle_field_solution;

  /// Finite element for the void fraction
  std::shared_ptr<FESystem<dim>> fe;

  /// Index set for the locally owned degree of freedoms
  IndexSet locally_owned_dofs;

  /// Index set for the locally relevant degree of freedoms
  IndexSet locally_relevant_dofs;

  /// Locally owned solution of the void fraction
  GlobalVectorType particle_field_locally_owned;

  /// System matrix used to assemble the smoothed L2 projection of the void
  /// fraction
  TrilinosWrappers::SparseMatrix system_matrix;

  /// Right-hand side used to assemble the smoothed L2 projection of the void
  /// fraction
  GlobalVectorType system_rhs;

  /// Constraints used for the boundary conditions of the void fraction.
  /// Currently, this is only used to establish periodic void fractions. This
  /// object has to be made public because the boundary conditions are set
  /// outside of the object for now.
  AffineConstraints<double> particle_field_constraints;

  /**
   * @brief Setup the degrees of freedom.
   *
   */
  void
  setup_dofs();
};

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
    , particle_handler(particle_handler)
    , particle_velocity_qcm(triangulation, fe_degree, simplex)
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
        if (this->void_fraction_parameters->quadrature_rule ==
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
      }
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
   * @brief Assemble and solve the system.
   */
  void
  solve() override
  {}


  /**
   * @brief Percolate the time vector for the void fraction. This operation is called at the end of a time-step.
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

  /// deal.II vector for the void fraction
  LinearAlgebra::distributed::Vector<double> void_fraction_solution;

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

  /*
   *
   */
  template <int component_start, int n_components>
  void
  calculate_field_projection(
    ParticleFieldQCM<dim, component_start, n_components> &field_qcm)
  {
    FEValues<dim> fe_values_velocity(*mapping,
                                     *field_qcm.fe,
                                     *quadrature,
                                     update_values | update_quadrature_points |
                                       update_JxW_values | update_gradients);

    FEValuesExtractors::Vector velocities;
    velocities.first_vector_component = 0;

    const unsigned int dofs_per_cell = field_qcm.fe->dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const unsigned int                   n_q_points = quadrature->size();
    FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>              local_rhs(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_vf(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_vf(dofs_per_cell);

    double         r_sphere = 0.0;
    double         total_volume_of_particle_in_sphere;
    Tensor<1, dim> particles_velocity_in_sphere;
    double qcm_sphere_diameter = void_fraction_parameters->qcm_sphere_diameter;

    // If the reference sphere diameter is user-defined, the radius is
    // calculated from it, otherwise, the value must be calculated while looping
    // over the cells.
    bool calculate_reference_sphere_radius = true;
    if (qcm_sphere_diameter > 1e-16)
      {
        r_sphere                          = 0.5 * qcm_sphere_diameter;
        calculate_reference_sphere_radius = false;
      }

    // Lambda functions for calculating the radius of the reference sphere
    // Calculate the radius by the volume (area in 2D) of sphere:
    // r = (2*dim*V/pi)^(1/dim) / 2
    auto radius_sphere_volume_cell = [](auto cell_measure) {
      return 0.5 * pow(2.0 * dim * cell_measure / M_PI, 1.0 / double(dim));
    };

    // Calculate the radius is obtained from the volume of sphere based on
    // R_s = h_omega:
    // V_s = pi*(2*V_c^(1/dim))^(dim)/(2*dim) = pi*2^(dim)*V_c/(2*dim)
    auto radius_h_omega = [&radius_sphere_volume_cell](double cell_measure) {
      double reference_sphere_volume =
        M_PI * Utilities::fixed_power<dim>(2.0) * cell_measure / (2.0 * dim);

      return radius_sphere_volume_cell(reference_sphere_volume);
    };

    field_qcm.system_rhs    = 0;
    field_qcm.system_matrix = 0;


    for (const auto &cell : field_qcm.dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values_velocity.reinit(cell);

            local_matrix = 0;
            local_rhs    = 0;

            // Array of real locations for the quadrature points
            std::vector<Point<dim>> quadrature_point_location;

            quadrature_point_location =
              fe_values_velocity.get_quadrature_points();

            // Active neighbors include the current cell as well
            auto active_neighbors =
              LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell,
                                                          cell);

            // Periodic neighbors of the current cell
            auto active_periodic_neighbors =
              LetheGridTools::find_cells_around_cell<dim>(
                vertices_to_periodic_cell, cell);

            // Define the volume of the reference sphere to be used as the
            // averaging volume for the QCM
            if (calculate_reference_sphere_radius)
              {
                if (void_fraction_parameters->qcm_sphere_equal_cell_volume ==
                    true)
                  {
                    // Get the radius by the volume of sphere which is
                    // equal to the volume of cell
                    r_sphere = radius_sphere_volume_cell(cell->measure());
                  }
                else
                  {
                    // The radius is obtained from the volume of sphere based
                    // on R_s = h_omega
                    r_sphere = radius_h_omega(cell->measure());
                  }
              }

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                total_volume_of_particle_in_sphere = 0;
                particles_velocity_in_sphere       = 0;

                for (unsigned int m = 0; m < active_neighbors.size(); m++)
                  {
                    // Loop over particles in neighbor cell
                    // Begin and end iterator for particles in neighbor cell
                    const auto pic =
                      particle_handler->particles_in_cell(active_neighbors[m]);
                    for (auto &particle : pic)
                      {
                        double distance            = 0;
                        auto   particle_properties = particle.get_properties();
                        const double r_particle =
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                          0.5;
                        double single_particle_volume =
                          M_PI * Utilities::fixed_power<dim>(r_particle * 2.0) /
                          (2 * dim);

                        // Distance between particle and quadrature point
                        // centers
                        distance = particle.get_location().distance(
                          quadrature_point_location[q]);

                        // Particle completely in the reference sphere
                        if (distance <= (r_sphere - r_particle))
                          {
                            const double particle_volume =
                              single_particle_volume;

                            total_volume_of_particle_in_sphere +=
                              particle_volume;
                            for (unsigned int d = 0; d < n_components; ++d)
                              {
                                particles_velocity_in_sphere[d] +=
                                  particle_volume *
                                  particle_properties[component_start + d];
                              }
                          }

                        // Particle partially in the reference sphere
                        else if ((distance > (r_sphere - r_particle)) &&
                                 (distance < (r_sphere + r_particle)))
                          {
                            if (dim == 2)
                              {
                                const double particle_volume =
                                  particle_circle_intersection_2d(r_particle,
                                                                  r_sphere,
                                                                  distance);

                                total_volume_of_particle_in_sphere +=
                                  particle_volume;

                                for (unsigned int d = 0; d < n_components; ++d)
                                  {
                                    particles_velocity_in_sphere[d] +=
                                      particle_volume *
                                      particle_properties[component_start + d];
                                  }
                              }
                            else if (dim == 3)
                              {
                                const double particle_volume =
                                  particle_sphere_intersection_3d(r_particle,
                                                                  r_sphere,
                                                                  distance);
                                total_volume_of_particle_in_sphere +=
                                  particle_volume;

                                for (unsigned int d = 0; d < n_components; ++d)
                                  {
                                    particles_velocity_in_sphere[d] +=
                                      particle_volume *
                                      particle_properties[component_start + d];
                                  }
                              }
                          }

                        // Particle completely outside the reference sphere. Do
                        // absolutely nothing.
                      }
                  }

                // Execute same operations for periodic neighbors, if the
                // simulation has no periodic boundaries, the container is
                // empty. Also, those operations cannot be done in the previous
                // loop because the particles on the periodic side need a
                // correction with an offset for the distance with the
                // quadrature point
                for (unsigned int m = 0; m < active_periodic_neighbors.size();
                     m++)
                  {
                    // Loop over particles in periodic neighbor cell
                    const auto pic = particle_handler->particles_in_cell(
                      active_periodic_neighbors[m]);
                    for (auto &particle : pic)
                      {
                        double distance            = 0;
                        auto   particle_properties = particle.get_properties();
                        const double r_particle =
                          particle_properties
                            [DEM::CFDDEMProperties::PropertiesIndex::dp] *
                          0.5;
                        double single_particle_volume =
                          M_PI * Utilities::fixed_power<dim>(r_particle * 2) /
                          (2 * dim);

                        // Adjust the location of the particle in the cell to
                        // account for the periodicity. If the position of the
                        // periodic cell if greater than the position of the
                        // current cell, the particle location needs a negative
                        // correction, and vice versa. Since the particle is in
                        // the periodic cell, this correction is the inverse of
                        // the correction for the volumetric contribution
                        const Point<dim> particle_location =
                          (active_periodic_neighbors[m]
                             ->center()[periodic_direction] >
                           cell->center()[periodic_direction]) ?
                            particle.get_location() - periodic_offset :
                            particle.get_location() + periodic_offset;

                        // Distance between particle and quadrature point
                        // centers
                        distance = particle_location.distance(
                          quadrature_point_location[q]);

                        // Particle completely in the reference sphere
                        if (distance <= (r_sphere - r_particle))
                          {
                            const double particle_volume =
                              single_particle_volume;
                            total_volume_of_particle_in_sphere +=
                              particle_volume;
                            for (unsigned int d = 0; d < n_components; ++d)
                              {
                                particles_velocity_in_sphere[d] +=
                                  particle_volume *
                                  particle_properties[component_start + d];
                              }
                          }

                        // Particle partially in the reference sphere
                        else if ((distance > (r_sphere - r_particle)) &&
                                 (distance < (r_sphere + r_particle)))
                          {
                            if (dim == 2)
                              {
                                const double particle_volume =
                                  particle_circle_intersection_2d(r_particle,
                                                                  r_sphere,
                                                                  distance);

                                total_volume_of_particle_in_sphere +=
                                  particle_volume;

                                for (unsigned int d = 0; d < n_components; ++d)
                                  {
                                    particles_velocity_in_sphere[d] +=
                                      particle_volume *
                                      particle_properties[component_start + d];
                                  }
                              }
                            else if (dim == 3)
                              {
                                const double particle_volume =
                                  particle_sphere_intersection_3d(r_particle,
                                                                  r_sphere,
                                                                  distance);
                                total_volume_of_particle_in_sphere +=
                                  particle_volume;

                                for (unsigned int d = 0; d < dim; ++d)
                                  {
                                    particles_velocity_in_sphere[d] +=
                                      particle_volume *
                                      particle_properties
                                        [DEM::CFDDEMProperties::
                                           PropertiesIndex::v_x +
                                         d];
                                  }
                              }
                          }


                        // Particle completely outside the reference sphere. Do
                        // absolutely nothing.
                      }
                  }

                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                  {
                    phi_vf[k] = fe_values_velocity[velocities].value(k, q);
                    grad_phi_vf[k] =
                      fe_values_velocity[velocities].gradient(k, q);
                  }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    // Assemble L2 projection
                    // Matrix assembly
                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        local_matrix(i, j) +=
                          ((phi_vf[j] * phi_vf[i]) +
                           (this->l2_smoothing_factor *
                            scalar_product(grad_phi_vf[j], grad_phi_vf[i]))) *
                          fe_values_velocity.JxW(q);
                      }

                    if (total_volume_of_particle_in_sphere > 0)
                      {
                        local_rhs(i) += phi_vf[i] *
                                        particles_velocity_in_sphere /
                                        total_volume_of_particle_in_sphere *
                                        fe_values_velocity.JxW(q);
                      }
                  }
              }

            cell->get_dof_indices(local_dof_indices);
            field_qcm.field_constraints.distribute_local_to_global(
              local_matrix,
              local_rhs,
              local_dof_indices,
              field_qcm.system_matrix,
              field_qcm.system_rhs);
          }
      }

    field_qcm.system_matrix.compress(VectorOperation::add);
    field_qcm.system_rhs.compress(VectorOperation::add);

    // Solve the L2 projection system
    const double linear_solver_tolerance =
      linear_solver_parameters.minimum_residual;

    if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Tolerance of iterative solver is : "
                    << linear_solver_tolerance << std::endl;
      }

    SolverControl solver_control(linear_solver_parameters.max_iterations,
                                 linear_solver_tolerance,
                                 true,
                                 true);

    TrilinosWrappers::SolverCG solver(solver_control);

    //**********************************************
    // Trillinos Wrapper ILU Preconditioner
    //*********************************************
    const double ilu_fill = linear_solver_parameters.ilu_precond_fill;
    const double ilu_atol = linear_solver_parameters.ilu_precond_atol;
    const double ilu_rtol = linear_solver_parameters.ilu_precond_rtol;

    TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
      ilu_fill, ilu_atol, ilu_rtol, 0);

    ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

    ilu_preconditioner->initialize(field_qcm.system_matrix,
                                   preconditionerOptions);

    solver.solve(field_qcm.system_matrix,
                 field_qcm.field_locally_owned,
                 field_qcm.system_rhs,
                 *ilu_preconditioner);

    if (linear_solver_parameters.verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step() << " steps " << std::endl;
      }

    field_qcm.field_constraints.distribute(field_qcm.field_locally_owned);
    field_qcm.field_locally_relevant = field_qcm.field_locally_owned;

#ifndef LETHE_USE_LDV
    // Perform copy between two vector types to ensure there is a deal.II vector
    convert_vector_trilinos_to_dealii(field_qcm.field_solution,
                                      field_qcm.field_locally_relevant);
    field_qcm.field_solution.update_ghost_values();
#else
    void_fraction_solution = void_fraction_locally_relevant;
#endif
  }

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

private:
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

  // Smoothing length factor for the void fraction calculation
  const double l2_smoothing_factor =
    Utilities::fixed_power<2>(void_fraction_parameters->l2_smoothing_length);

public:
  ParticleFieldQCM<dim, DEM::CFDDEMProperties::PropertiesIndex::v_x, 3>
    particle_velocity_qcm;
};



#endif
