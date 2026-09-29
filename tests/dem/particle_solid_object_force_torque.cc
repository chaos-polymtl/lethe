// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Check the force and torque exerted by the particles on solid surfaces.
 *
 * Three moving and rotating particles are in contact with two flat solid
 * surfaces: two particles on opposite sides of the first surface, and one
 * particle on the second surface. The contacts are off-center, and the
 * friction and rolling resistance are active, so that the normal force, the
 * tangential force, the tangential torque and the rolling resistance torque
 * all contribute. For every solid, the test verifies that:
 * - the force on the solid is the opposite of the contact force applied on
 * the particles in contact with it (action-reaction);
 * - the torque on the solid about its center of rotation is the opposite of
 * the moment, about the same point, of the force and torque applied on these
 * particles (angular momentum balance).
 * The contact forces are evaluated twice to verify that the loads on the
 * solids are those of the last evaluation and are not accumulated.
 */

// Deal.II
#include <deal.II/base/function.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters.h>
#include <core/serial_solid.h>
#include <core/solid_objects_parameters.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>

// Tests (with common definitions)
#include <../tests/dem/test_particles_functions.h>

#include <../tests/tests.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Create a square solid surface made of triangles, initially lying in
 * the z = 0 plane and centered at the origin, then rotated and translated.
 *
 * @param[in] half_width Half of the width of the square.
 * @param[in] rotation_axis Axis of the initial rotation of the surface.
 * @param[in] rotation_angle Angle of the initial rotation of the surface.
 * @param[in] translation Initial translation of the surface.
 * @param[in] center_of_rotation Center of rotation of the solid.
 * @param[in] id Identifier of the solid.
 *
 * @return The solid surface.
 */
std::shared_ptr<SerialSolid<2, 3>>
create_square_solid_surface(const double        half_width,
                            const Tensor<1, 3> &rotation_axis,
                            const double        rotation_angle,
                            const Tensor<1, 3> &translation,
                            const Point<3>     &center_of_rotation,
                            const unsigned int  id)
{
  auto param = std::make_shared<Parameters::RigidSolidObject<2, 3>>();

  param->solid_mesh.type           = Parameters::Mesh<2, 3>::Type::dealii;
  param->solid_mesh.grid_type      = "hyper_cube";
  param->solid_mesh.grid_arguments = std::to_string(-half_width) + " : " +
                                     std::to_string(half_width) + " : false";
  param->solid_mesh.initial_refinement = 0;
  param->solid_mesh.simplex            = true;
  param->solid_mesh.rotation_axis      = rotation_axis;
  param->solid_mesh.rotation_angle     = rotation_angle;
  param->solid_mesh.translation        = translation;
  param->center_of_rotation            = center_of_rotation;
  param->output_bool                   = false;
  param->thermal_boundary_type = Parameters::ThermalBoundaryType::adiabatic;
  param->translational_velocity =
    std::make_shared<Functions::ZeroFunction<3>>(3);
  param->angular_velocity  = std::make_shared<Functions::ZeroFunction<3>>(3);
  param->solid_temperature = std::make_shared<Functions::ZeroFunction<3>>(1);

  return std::make_shared<SerialSolid<2, 3>>(param, id);
}

/**
 * @brief Verify the action-reaction and the angular momentum balance between
 * the particles and each solid, and print the loads on the solids.
 *
 * @param[in] particle_handler Particle handler.
 * @param[in] contact_outcome Force and torque applied on the particles.
 * @param[in] solids Solid surfaces.
 * @param[in] particle_solid Index of the solid in contact with each particle,
 * indexed by particle id.
 * @param[in] solid_forces Force exerted by the particles on each solid.
 * @param[in] solid_torques Torque exerted by the particles on each solid about
 * its center of rotation.
 */
template <typename PropertiesIndex>
void
check_balance(
  const Particles::ParticleHandler<3>                   &particle_handler,
  const ParticleInteractionOutcomes<PropertiesIndex>    &contact_outcome,
  const std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solids,
  const std::vector<unsigned int>                       &particle_solid,
  const std::vector<Tensor<1, 3>>                       &solid_forces,
  const std::vector<Tensor<1, 3>>                       &solid_torques)
{
  // Relative tolerance on the balances. They are satisfied up to round-off
  // since the same contact forces are applied on the particles and on the
  // solids.
  const double balance_tolerance = 1e-12;

  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      const Point<3> center_of_rotation =
        solids[i_solid]->get_center_of_rotation();

      // Sum of the loads on the solid and on the particles in contact with it.
      // The moment of the loads on the particles is taken about the center of
      // rotation of the solid.
      Tensor<1, 3> force_balance  = solid_forces[i_solid];
      Tensor<1, 3> torque_balance = solid_torques[i_solid];
      Tensor<1, 3> particle_direct_torque;
      for (const auto &particle : particle_handler)
        {
          if (particle_solid[particle.get_id()] != i_solid)
            continue;

          const unsigned int  local_index = particle.get_local_index();
          const Tensor<1, 3> &force       = contact_outcome.force[local_index];
          const Tensor<1, 3> &torque      = contact_outcome.torque[local_index];

          force_balance += force;
          torque_balance +=
            cross_product_3d(particle.get_location() - center_of_rotation,
                             force) +
            torque;
          particle_direct_torque += torque;
        }

      deallog << "Solid " << i_solid << std::endl;
      deallog << "  Force on solid: " << solid_forces[i_solid] << std::endl;
      deallog << "  Torque on solid: " << solid_torques[i_solid] << std::endl;
      deallog << "  Tangential and rolling resistance torques on particles: "
              << particle_direct_torque << std::endl;
      deallog << "  Action-reaction satisfied: "
              << (force_balance.norm() <=
                      balance_tolerance * solid_forces[i_solid].norm() ?
                    "true" :
                    "false")
              << std::endl;
      deallog << "  Angular momentum balance satisfied: "
              << (torque_balance.norm() <=
                      balance_tolerance * solid_torques[i_solid].norm() ?
                    "true" :
                    "false")
              << std::endl;
    }
}

template <typename PropertiesIndex,
          Parameters::Lagrangian::ParticleWallContactForceModel contact_model,
          Parameters::Lagrangian::RollingResistanceMethod       rolling_model>
void
test()
{
  // Background triangulation holding the particles
  parallel::distributed::Triangulation<3> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  triangulation.refine_global(2);
  MappingQ<3> mapping(1);

  // DEM parameters
  DEMSolverParameters<3> dem_parameters;
  set_default_dem_parameters(1, dem_parameters);
  auto &properties = dem_parameters.lagrangian_physical_properties;
  properties.youngs_modulus_particle[0]                      = 1e7;
  properties.youngs_modulus_wall                             = 1e7;
  properties.poisson_ratio_particle[0]                       = 0.3;
  properties.poisson_ratio_wall                              = 0.3;
  properties.restitution_coefficient_particle[0]             = 0.5;
  properties.restitution_coefficient_wall                    = 0.5;
  properties.friction_coefficient_particle[0]                = 0.4;
  properties.friction_coefficient_wall                       = 0.4;
  properties.rolling_friction_coefficient_particle[0]        = 0.2;
  properties.rolling_friction_wall                           = 0.2;
  properties.rolling_viscous_damping_coefficient_particle[0] = 0.2;
  properties.rolling_viscous_damping_wall                    = 0.2;
  dem_parameters.model_parameters.rolling_resistance_method  = rolling_model;
  const double dt                                            = 1e-5;
  const double particle_diameter                             = 0.01;
  const double particle_mass                                 = 1e-3;

  // Solid 0 is a horizontal square in the z = 0 plane. Solid 1 is a vertical
  // square in the x = 0.7 plane. Their centers of rotation are away from the
  // contacts so that the lever arms are not trivial.
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> solids;
  solids.push_back(create_square_solid_surface(0.5,
                                               Tensor<1, 3>({1., 0., 0.}),
                                               0.,
                                               Tensor<1, 3>(),
                                               Point<3>(0.3, -0.2, 0.4),
                                               0));
  solids.push_back(create_square_solid_surface(0.2,
                                               Tensor<1, 3>({0., 1., 0.}),
                                               0.5 * numbers::PI,
                                               Tensor<1, 3>({0.7, 0., 0.}),
                                               Point<3>(0.9, 0.3, -0.1),
                                               1));

  // Particles 0 and 1 overlap solid 0 from above and from below, and particle
  // 2 overlaps solid 1. The overlaps are a few percent of the diameter.
  std::vector<Point<3>>     positions{Point<3>(0.2, 0.1, 0.0046),
                                  Point<3>(-0.15, -0.3, -0.0047),
                                  Point<3>(0.6955, 0.05, 0.03)};
  std::vector<Tensor<1, 3>> velocities{Tensor<1, 3>({0.3, -0.1, -0.2}),
                                       Tensor<1, 3>({-0.2, 0.05, 0.15}),
                                       Tensor<1, 3>({0.25, 0.1, -0.1})};
  std::vector<Tensor<1, 3>> angular_velocities{Tensor<1, 3>({4., -2., 3.}),
                                               Tensor<1, 3>({-3., 1., 2.}),
                                               Tensor<1, 3>({1., 5., -2.})};
  const std::vector<unsigned int> particle_solid{0, 0, 1};

  Particles::ParticleHandler<3> particle_handler(triangulation,
                                                 mapping,
                                                 PropertiesIndex::n_properties);
  for (unsigned int i = 0; i < positions.size(); ++i)
    {
      Particles::ParticleIterator<3> particle = construct_particle_iterator<3>(
        particle_handler, triangulation, positions[i], i);
      set_particle_properties<3, PropertiesIndex>(particle,
                                                  0,
                                                  particle_diameter,
                                                  particle_mass,
                                                  velocities[i],
                                                  angular_velocities[i]);
    }
  particle_handler.sort_particles_into_subdomains_and_cells();

  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  contact_outcome.resize_interaction_containers(
    particle_handler.get_max_local_particle_index());

  // Every particle is a contact candidate of every triangle of every solid.
  // The contact force calculation keeps only the actual contacts.
  typename DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    contact_candidates(solids.size());
  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    for (const auto &triangle :
         solids[i_solid]->get_triangulation()->active_cell_iterators())
      for (auto particle = particle_handler.begin();
           particle != particle_handler.end();
           ++particle)
        contact_candidates[i_solid][triangle].emplace(particle->get_id(),
                                                      particle);

  typename DEM::dem_data_structures<
    3>::particle_floating_mesh_potentially_in_contact potentially_in_contact;
  particle_floating_mesh_fine_search<3>(contact_candidates,
                                        potentially_in_contact);

  ParticleWallContactForce<3, PropertiesIndex, contact_model, rolling_model>
    force_object(dem_parameters);

  std::vector<Tensor<1, 3>> solid_forces;
  std::vector<Tensor<1, 3>> solid_torques;

  for (unsigned int evaluation = 0; evaluation < 2; ++evaluation)
    {
      reinitialize_contact_outcomes<3, PropertiesIndex>(particle_handler,
                                                        contact_outcome);

      force_object.calculate_particle_solid_object_contact(
        potentially_in_contact,
        dt,
        solids,
        contact_outcome,
        solid_forces,
        solid_torques);

      deallog << "Contact force evaluation " << evaluation << std::endl;
      check_balance(particle_handler,
                    contact_outcome,
                    solids,
                    particle_solid,
                    solid_forces,
                    solid_torques);
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();

      using namespace Parameters::Lagrangian;

      deallog << "Non-linear contact model and constant rolling resistance"
              << std::endl;
      test<DEM::DEMProperties::PropertiesIndex,
           ParticleWallContactForceModel::nonlinear,
           RollingResistanceMethod::constant>();

      deallog << "Linear contact model and viscous rolling resistance"
              << std::endl;
      test<DEM::DEMProperties::PropertiesIndex,
           ParticleWallContactForceModel::linear,
           RollingResistanceMethod::viscous>();
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
