// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the tangential force of the contact model of the sharp immersed
 * boundary DEM between two bodies made of different materials. Without normal
 * relative velocity and without contact history, the tangential force must be
 * 8 G* sqrt(R* delta) v_t dt, where G* is the effective shear modulus used by
 * the DEM, and it must not depend on the order of the two bodies.
 */

// Deal.II includes
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/tests.h>

#include <cmath>

using namespace dealii;

template <int dim>
void
test()
{
  using ObjectProperties = typename IBParticlesDEM<dim>::ObjectProperties;
  using ContactInfo      = typename IBParticlesDEM<dim>::ContactInfo;

  ObjectProperties body_one;
  body_one.object_mass                         = 1.;
  body_one.object_radius                       = 0.1;
  body_one.object_youngs_modulus               = 1e7;
  body_one.object_poisson_ratio                = 0.2;
  body_one.object_restitution_coefficient      = 0.9;
  body_one.object_friction_coefficient         = 0.5;
  body_one.object_rolling_friction_coefficient = 0.;

  ObjectProperties body_two;
  body_two.object_mass                         = 2.;
  body_two.object_radius                       = 0.2;
  body_two.object_youngs_modulus               = 3e7;
  body_two.object_poisson_ratio                = 0.35;
  body_two.object_restitution_coefficient      = 0.9;
  body_two.object_friction_coefficient         = 0.5;
  body_two.object_rolling_friction_coefficient = 0.;

  const double normal_overlap      = 1e-3;
  const double tangential_velocity = 0.1;
  const double dt                  = 1e-5;
  const double tolerance           = 1e-12;

  // Tangential force when the first body slides along x on the second one.
  // The normal of the contact is along y, so the relative velocity is purely
  // tangential. The force is small enough to stay below the Coulomb limit.
  auto compute_tangential_force = [&](const ObjectProperties &first,
                                      const ObjectProperties &second) {
    IBParticlesDEM<dim> ib_dem;
    ContactInfo         contact_info{};
    Point<3>            contact_point;
    Tensor<1, 3>        contact_normal;
    contact_normal[1] = 1.;
    Tensor<1, 3> normal_force;
    Tensor<1, 3> tangential_force;
    Tensor<1, 3> rolling_resistance_torque;
    Point<dim>   position_one;
    position_one[1] = -first.object_radius;
    Point<dim> position_two;
    position_two[1] = second.object_radius;
    Tensor<1, 3> velocity_one;
    velocity_one[0] = tangential_velocity;
    Tensor<1, 3> velocity_two;
    Tensor<1, 3> omega_one;
    Tensor<1, 3> omega_two;

    ib_dem.calculate_force_model(normal_overlap,
                                 contact_info,
                                 contact_point,
                                 contact_normal,
                                 normal_force,
                                 tangential_force,
                                 rolling_resistance_torque,
                                 position_one,
                                 velocity_one,
                                 omega_one,
                                 first,
                                 position_two,
                                 velocity_two,
                                 omega_two,
                                 second,
                                 dt);
    return tangential_force;
  };

  // Effective properties of the DEM (see ParticleWallContactForce).
  const double E_1  = body_one.object_youngs_modulus;
  const double E_2  = body_two.object_youngs_modulus;
  const double nu_1 = body_one.object_poisson_ratio;
  const double nu_2 = body_two.object_poisson_ratio;
  const double effective_shear_modulus =
    E_1 * E_2 /
    (2. * (E_2 * (2. - nu_1) * (1. + nu_1) + E_1 * (2. - nu_2) * (1. + nu_2)));
  const double effective_radius =
    body_one.object_radius * body_two.object_radius /
    (body_one.object_radius + body_two.object_radius);

  Tensor<1, 3> expected_tangential_force;
  expected_tangential_force[0] = 8. * effective_shear_modulus *
                                 std::sqrt(effective_radius * normal_overlap) *
                                 tangential_velocity * dt;

  const Tensor<1, 3> force_one_on_two =
    compute_tangential_force(body_one, body_two);
  const Tensor<1, 3> force_two_on_one =
    compute_tangential_force(body_two, body_one);

  const bool force_is_right =
    (force_one_on_two - expected_tangential_force).norm() <
    tolerance * expected_tangential_force.norm();
  const bool force_is_symmetric = (force_one_on_two - force_two_on_one).norm() <
                                  tolerance * expected_tangential_force.norm();

  deallog << "Dimension " << dim << std::endl;
  deallog << "Tangential force matches the effective shear modulus of the DEM: "
          << (force_is_right ? "true" : "false") << std::endl;
  deallog << "Tangential force does not depend on the order of the bodies: "
          << (force_is_symmetric ? "true" : "false") << std::endl;
}

int
main()
{
  try
    {
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
