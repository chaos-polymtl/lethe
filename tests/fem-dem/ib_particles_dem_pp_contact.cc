// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the normal force of the particle-particle contact of the sharp
 * immersed boundary DEM. Two particles at rest overlap by a known amount
 * delta. Without velocity, the contact force must be the Hertz force of the
 * full overlap, 4/3 E* sqrt(R*) delta^(3/2), along the line joining the
 * particles and pushing them apart.
 *
 * Sphere-sphere: two spheres of different radii, whose line of centers is not
 * aligned with an axis.
 *
 * Sphere-box: a sphere overlaps a face of a box. The force on the sphere must
 * not depend on the order of the two particles, since the distance routine of
 * the sphere is used in both cases and the normal it returns is flipped when
 * the sphere is the second particle.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/ib_particle.h>
#include <core/parameters.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <cmath>
#include <memory>
#include <vector>

using namespace dealii;

constexpr double dt_dem         = 1e-5;
constexpr double normal_overlap = 1e-3;

// The DEM uses the constant 0.66665 in place of 2/3 in the normal spring
// constant, which differs from the Hertz force by 2.5e-5 in relative terms.
constexpr double hertz_tolerance = 1e-4;

/**
 * @brief Compute the particle-particle contact force on each particle.
 *
 * @param[in] parameters Parameters of the IB particles.
 * @param[in] particles Particles, at rest.
 *
 * @return The contact force on each particle.
 */
template <int dim>
std::vector<Tensor<1, 3>>
particle_particle_contact_force(
  const std::shared_ptr<Parameters::IBParticles<dim>> &parameters,
  const std::vector<IBParticle<dim>>                  &particles)
{
  IBParticlesDEM<dim> ib_dem;
  ib_dem.initialize(parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
  ib_dem.update_contact_candidates();

  std::vector<Tensor<1, 3>> force(particles.size());
  std::vector<Tensor<1, 3>> torque(particles.size());
  ib_dem.calculate_pp_contact_force(dt_dem, force, torque);
  return force;
}

/**
 * @brief Hertz normal force between two particles made of the same material,
 * for the overlap normal_overlap.
 *
 * @param[in] particle_one First particle.
 * @param[in] particle_two Second particle.
 *
 * @return The magnitude of the force.
 */
template <int dim>
double
hertz_force(const IBParticle<dim> &particle_one,
            const IBParticle<dim> &particle_two)
{
  const double effective_youngs_modulus =
    particle_one.youngs_modulus /
    (2. * (1. - particle_one.poisson_ratio * particle_one.poisson_ratio));
  const double effective_radius = particle_one.radius * particle_two.radius /
                                  (particle_one.radius + particle_two.radius);
  return 4. / 3. * effective_youngs_modulus * std::sqrt(effective_radius) *
         std::pow(normal_overlap, 1.5);
}

/**
 * @brief Build a box at rest whose motion is integrated, with the material of
 * the spheres of make_ib_sphere.
 *
 * @param[in] id Identifier of the particle.
 * @param[in] position Position of the center of the box.
 *
 * @return The box, whose half lengths are 0.2, 0.3 and 0.25.
 */
template <int dim>
IBParticle<dim>
make_ib_box(const unsigned int id, const Point<dim> &position)
{
  IBParticle<dim>     particle = make_ib_sphere<dim>(id, 0.1, position);
  std::vector<double> half_lengths{0.2, 0.3, 0.25};
  half_lengths.resize(dim);
  particle.initialize_shape("hyper rectangle", half_lengths);
  particle.radius = particle.shape->effective_radius;
  particle.set_position(position);
  particle.set_orientation(particle.orientation);
  particle.initialize_previous_solution();
  return particle;
}

template <int dim>
void
test()
{
  const double tolerance = 1e-12;
  auto         parameters =
    make_ib_particles_parameters<dim>((dim == 2) ? "0; 0" : "0; 0; 0");
  // Use the radius of the particles in the Hertz model rather than the
  // curvature of their surface at the contact point.
  parameters->approximate_radius_for_contact = true;

  auto to_string = [](const bool value) { return value ? "true" : "false"; };

  deallog << "Dimension " << dim << std::endl;

  // Sphere-sphere contact. The distance between the centers is
  // radius_one + radius_two - normal_overlap.
  {
    const double radius_one = 0.1;
    const double radius_two = 0.25;

    // Unit vector from the first sphere to the second one.
    Tensor<1, dim> direction;
    if constexpr (dim == 2)
      {
        direction[0] = 0.6;
        direction[1] = 0.8;
      }
    else
      {
        direction[0] = 2. / 3.;
        direction[1] = 1. / 3.;
        direction[2] = 2. / 3.;
      }

    Point<dim> position_one;
    for (int d = 0; d < dim; ++d)
      position_one[d] = 0.1 * (d + 1);
    const Point<dim> position_two =
      position_one + (radius_one + radius_two - normal_overlap) * direction;

    const std::vector<IBParticle<dim>> particles{
      make_ib_sphere<dim>(0, radius_one, position_one),
      make_ib_sphere<dim>(1, radius_two, position_two)};
    const auto force =
      particle_particle_contact_force<dim>(parameters, particles);

    const Tensor<1, 3> expected_force =
      -hertz_force(particles[0], particles[1]) * tensor_nd_to_3d(direction);

    deallog << "Sphere-sphere contact force is the Hertz force of the full "
               "overlap: "
            << to_string((force[0] - expected_force).norm() <
                         hertz_tolerance * expected_force.norm())
            << std::endl;
    deallog << "Sphere-sphere contact forces are opposite: "
            << to_string((force[0] + force[1]).norm() <
                         tolerance * expected_force.norm())
            << std::endl;
  }

  // Sphere-box contact. The box is centered at the origin and the sphere
  // overlaps its face x = 0.2.
  {
    const double radius = 0.1;
    Point<dim>   sphere_position;
    sphere_position[0] = 0.2 + radius - normal_overlap;

    const auto force_sphere_first = particle_particle_contact_force<dim>(
      parameters,
      {make_ib_sphere<dim>(0, radius, sphere_position),
       make_ib_box<dim>(1, Point<dim>())});
    const std::vector<IBParticle<dim>> box_first{
      make_ib_box<dim>(0, Point<dim>()),
      make_ib_sphere<dim>(1, radius, sphere_position)};
    const auto force_box_first =
      particle_particle_contact_force<dim>(parameters, box_first);

    // The sphere is pushed away from the box, along x.
    Tensor<1, 3> expected_force;
    expected_force[0] = hertz_force(box_first[1], box_first[0]);

    deallog << "Sphere-box contact force is the Hertz force of the full "
               "overlap: "
            << to_string((force_sphere_first[0] - expected_force).norm() <
                         hertz_tolerance * expected_force.norm())
            << std::endl;
    deallog << "Sphere-box contact force does not depend on the order of the "
               "particles: "
            << to_string((force_sphere_first[0] - force_box_first[1]).norm() <
                           tolerance * expected_force.norm() &&
                         (force_sphere_first[1] - force_box_first[0]).norm() <
                           tolerance * expected_force.norm())
            << std::endl;
  }
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
