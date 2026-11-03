// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ib_particles_dem_test_utilities_h
#define lethe_ib_particles_dem_test_utilities_h

/**
 * @brief Utilities shared by the unit tests of IBParticlesDEM, the DEM used by
 * the sharp immersed boundary solver. They build the spherical particles and
 * the parameters of the DEM without parsing a parameter file.
 */

// Deal.II includes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/ib_particle.h>
#include <core/parameters.h>

#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Build a spherical IB particle at rest whose motion is integrated.
 *
 * @param[in] id Identifier of the particle, which must be its index in the
 * vector of particles given to IBParticlesDEM.
 * @param[in] radius Radius of the sphere.
 * @param[in] position Position of the center of the sphere.
 * @param[in] density Density of the particle.
 *
 * @return The particle, with its previous states equal to its current state.
 */
template <int dim>
IBParticle<dim>
make_ib_sphere(const unsigned int id,
               const double       radius,
               const Point<dim>  &position,
               const double       density = 1.)
{
  IBParticle<dim> particle;
  particle.initialize_all();
  particle.particle_id = id;
  particle.initialize_shape("sphere", std::vector<double>{radius});
  particle.radius = particle.shape->effective_radius;

  // The volume of a particle is its area in 2D.
  particle.volume = (dim == 2) ?
                      std::numbers::pi * radius * radius :
                      4. / 3. * std::numbers::pi * radius * radius * radius;
  particle.mass   = density * particle.volume;

  particle.youngs_modulus               = 1e8;
  particle.poisson_ratio                = 0.3;
  particle.restitution_coefficient      = 0.9;
  particle.friction_coefficient         = 0.;
  particle.rolling_friction_coefficient = 0.;
  particle.surface_energy               = 0.;
  particle.integrate_motion             = true;
  particle.mesh_based_precalculations   = false;

  particle.set_position(position);
  particle.set_orientation(particle.orientation);
  particle.initialize_previous_solution();
  return particle;
}

/**
 * @brief Build the parameters of the IB particles used by IBParticlesDEM.
 *
 * @param[in] gravity Expression of the gravity function, one component per
 * dimension separated by semicolons.
 *
 * @return The parameters, with their default values except for the gravity
 * and the lubrication force, which is disabled.
 */
template <int dim>
std::shared_ptr<Parameters::IBParticles<dim>>
make_ib_particles_parameters(const std::string &gravity)
{
  auto parameters = std::make_shared<Parameters::IBParticles<dim>>();
  parameters->enable_lubrication_force = false;

  ParameterHandler prm;
  Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
  prm.set("Function expression", gravity);
  parameters->f_gravity->parse_parameters(prm);
  return parameters;
}

#endif
