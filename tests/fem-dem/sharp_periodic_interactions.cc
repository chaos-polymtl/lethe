// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Compare seam interactions with the same physical pair in an interior frame.
#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/shape.h>

#include <fem-dem/ib_particles_dem.h>
#include <fem-dem/sharp_ib_periodic_boundaries.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <../tests/tests.h>

#include <cmath>
#include <exception>
#include <fstream>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <tuple>
#include <vector>

using namespace dealii;

/**
 * @brief Check periodic geometry, configuration and interaction invariants.
 * @tparam dim Spatial dimension.
 */
template <int dim>
void
test()
{
  Parameters::Lagrangian::BCDEM          dem;
  BoundaryConditions::BoundaryConditions cfd;
  dem.periodic_boundaries[0] = {.neighbor_id = 1, .direction = 0};
  cfd.type[0]                = BoundaryConditions::BoundaryType::periodic;
  cfd.periodic_boundaries[0] = {.neighbor_id = 1, .direction = 0};
  SharpIBPeriodicBoundaries<dim> helper;
  std::ostringstream             sink;
  ConditionalOStream             quiet(sink, false);
  helper.initialize(dem, cfd, quiet);
  Tensor<1, dim> offset;
  offset[0] = 1.;
  helper.set_periodic_offset(offset);
  for (const double x : {-4.2, -0.2, 0.2, 1.2, 5.2})
    {
      Point<dim> p;
      p[0]                  = x;
      const double expected = x - std::floor(x);
      AssertThrow(std::abs(helper.wrap_point(p, 0., 1.)[0] - expected) < 1e-12,
                  ExcMessage("Multi-period normalization failed"));
    }

  // Every accepted wrap must translate the complete BDF history uniformly.
  std::vector<Point<dim>> history(3);
  history[0][0]               = 0.8;
  history[1][0]               = 0.5;
  history[2][0]               = 0.2;
  const auto original_history = history;
  helper.translate_history(history, Tensor<1, dim>());
  AssertThrow(history == original_history,
              ExcMessage("History changed without a wrap"));
  Tensor<1, dim> wrap;
  wrap[0] = -2.;
  helper.translate_history(history, wrap);
  for (unsigned int i = 0; i < history.size(); ++i)
    AssertThrow(std::abs(history[i][0] - original_history[i][0] + 2.) < 1e-12,
                ExcMessage("History entries use different frames"));
  Point<dim> constituent_position;
  constituent_position[0] = -0.3;
  std::map<unsigned int, std::shared_ptr<Shape<dim>>> constituents;
  constituents[0] =
    std::make_shared<Sphere<dim>>(0.05, constituent_position, Tensor<1, 3>());
  constituent_position[0] = -0.28;
  constituents[1] =
    std::make_shared<Sphere<dim>>(0.05, constituent_position, Tensor<1, 3>());
  using Operation = typename CompositeShape<dim>::BooleanOperation;
  std::map<unsigned int, std::tuple<Operation, unsigned int, unsigned int>>
    operations;
  operations[2] = {Operation::Union, 0, 1};
  CompositeShape<dim> offset_shape(constituents,
                                   operations,
                                   Point<dim>(),
                                   Tensor<1, 3>());
  AssertThrow(!offset_shape.contains_rbf(),
              ExcMessage("Analytic composite incorrectly classified as RBF"));
  // A single node suffices: this checks shape-type validation, not RBF values.
  {
    std::ofstream input("periodic_rejected_rbf.input");
    input << "support_radius basis_function node_x node_y node_z weight\n"
          << "1 1 0 0 0 1\n";
  }
  auto rbf = std::make_shared<RBFShape<dim>>("periodic_rejected_rbf.input",
                                             Point<dim>(),
                                             Tensor<1, 3>());
  auto composite = std::make_shared<CompositeShape<dim>>(
    std::vector<std::shared_ptr<Shape<dim>>>{rbf},
    Point<dim>(),
    Tensor<1, 3>());
  CompositeShape<dim> nested(
    std::vector<std::shared_ptr<Shape<dim>>>{composite},
    Point<dim>(),
    Tensor<1, 3>());
  AssertThrow(composite->contains_rbf() && nested.contains_rbf(),
              ExcMessage("Composite RBF constituent was not detected"));
  Point<dim> center;
  center[0] = 0.3;
  AssertThrow(!helper
                 .generate_periodic_image_points(
                   center, helper.enclosing_radius(offset_shape), 0., 1.)
                 .empty(),
              ExcMessage("Offset geometry image was missed"));

  const auto distant_images =
    helper.generate_periodic_image_points(center, 2.4, 0., 1.);
  bool has_two_period_image = false;
  for (const auto &image : distant_images)
    has_two_period_image |= std::abs(image[0] - center[0] - 2.) < 1e-12;
  AssertThrow(has_two_period_image,
              ExcMessage("Multi-period geometric support was truncated"));

  // Invalid boundary pairs and zero-length periods must fail before use.
  bool rejected                          = false;
  cfd.periodic_boundaries[0].neighbor_id = 2;
  try
    {
      helper.initialize(dem, cfd, quiet);
    }
  catch (const std::exception &)
    {
      rejected = true;
    }
  AssertThrow(rejected,
              ExcMessage("Mismatched periodic configuration was accepted"));
  cfd.periodic_boundaries[0].neighbor_id = 1;
  helper.initialize(dem, cfd, quiet);
  rejected = false;
  try
    {
      helper.set_periodic_offset(Tensor<1, dim>());
    }
  catch (const std::exception &)
    {
      rejected = true;
    }
  AssertThrow(rejected, ExcMessage("Zero period was accepted"));

  typename IBParticlesDEM<dim>::ContactInfo sliding;
  sliding.normal_vector[0]              = 1.;
  sliding.tangential_displacement[0]    = 2.;
  sliding.tangential_displacement[1]    = 3.;
  sliding.tangential_displacement_limit = 0.25;
  sliding.limit_tangential_history();
  AssertThrow(std::abs(sliding.tangential_displacement[0]) < 1e-15 &&
                std::abs(sliding.tangential_displacement.norm() - 0.25) < 1e-15,
              ExcMessage(
                "Stored spring violates contact tangent or Coulomb limit"));

  auto parameters = std::make_shared<Parameters::IBParticles<dim>>();
  parameters->approximate_radius_for_contact = true;
  std::vector<IBParticle<dim>> particles(2);
  for (unsigned int i = 0; i < 2; ++i)
    {
      auto &p = particles[i];
      p.initialize_all();
      p.particle_id      = i;
      p.integrate_motion = true;
      p.radius           = 0.03;
      p.shape =
        std::make_shared<Sphere<dim>>(p.radius, p.position, p.orientation);
      p.mass                         = 1.;
      p.youngs_modulus               = 1e5;
      p.poisson_ratio                = 0.25;
      p.restitution_coefficient      = 0.8;
      p.friction_coefficient         = 0.3;
      p.rolling_friction_coefficient = 0.;
      p.position[0]                  = i == 0 ? 0.02 : 0.98;
      p.velocity[1]                  = i == 0 ? 0.01 : -0.01;
      p.set_position(p.position);
    }
  IBParticlesDEM<dim> seam, interior;
  auto walls = std::make_shared<Parameters::Lagrangian::FloatingWalls<dim>>();
  auto boundaries = std::make_shared<Parameters::Lagrangian::BCDEM>(dem);
  seam.initialize(parameters, walls, boundaries, MPI_COMM_WORLD, particles);
  seam.set_periodic_offset(offset);
  // Use independent shapes: changing the control must not mutate the seam pair.
  for (auto &p : particles)
    p.shape = p.shape->static_copy();
  particles[1].position[0] -= 1.;
  particles[1].set_position(particles[1].position);
  interior.initialize(parameters, walls, boundaries, MPI_COMM_WORLD, particles);
  // A distant seam provides the same periodic contact model in an interior
  // frame.
  interior.set_periodic_offset(10. * offset);
  seam.update_contact_candidates();
  interior.update_contact_candidates();
  for (unsigned int step = 0; step < 2; ++step)
    {
      // Changing the represented image must not reset physical pair history.
      if (step == 1)
        seam.dem_particles[0].set_position(1.02, 0);
      std::vector<Tensor<1, 3>> seam_force(2), interior_force(2),
        seam_torque(2), interior_torque(2);
      seam.calculate_pp_contact_force(1e-4, seam_force, seam_torque);
      interior.calculate_pp_contact_force(1e-4,
                                          interior_force,
                                          interior_torque);
      AssertThrow(seam_force[0].norm() > 0.,
                  ExcMessage("Seam contact was missed"));
      for (unsigned int i = 0; i < 2; ++i)
        {
          AssertThrow((seam_force[i] - interior_force[i]).norm() < 1e-9,
                      ExcMessage("Contact force differs"));
          AssertThrow((seam_torque[i] - interior_torque[i]).norm() < 1e-9,
                      ExcMessage("Contact torque differs"));
        }
      AssertThrow((seam_force[0] + seam_force[1]).norm() < 1e-12,
                  ExcMessage("Contact violates action/reaction"));
      AssertThrow(seam.dem_particles[1].position[0] > 0.9,
                  ExcMessage("Image evaluation changed primary geometry"));
    }
  AssertThrow(seam.get_particle_particle_contact_map()
                  .at(0)
                  .at(1)
                  .tangential_displacement.norm() > 0.,
              ExcMessage("Tangential history was not retained"));
  // Restart must preserve the accepted tangential spring.
  std::stringstream checkpoint;
  seam.write_contact_history(checkpoint);
  IBParticlesDEM<dim> restored;
  restored.read_contact_history(checkpoint);
  AssertThrow((restored.get_particle_particle_contact_map()
                 .at(0)
                 .at(1)
                 .tangential_displacement -
               seam.get_particle_particle_contact_map()
                 .at(0)
                 .at(1)
                 .tangential_displacement)
                  .norm() < 1e-15,
              ExcMessage("Contact checkpoint lost spring history"));
  // Repeating a nonlinear trial at the same time must not accumulate history.
  seam.update_particles(seam.dem_particles, 1., true);
  const auto accepted_spring = seam.get_particle_particle_contact_map()
                                 .at(0)
                                 .at(1)
                                 .tangential_displacement;
  std::vector<Tensor<1, 3>> trial_force(2), trial_torque(2);
  seam.calculate_pp_contact_force(1e-4, trial_force, trial_torque);
  const auto trial_spring = seam.get_particle_particle_contact_map()
                              .at(0)
                              .at(1)
                              .tangential_displacement;
  AssertThrow((trial_spring - accepted_spring).norm() > 1e-12,
              ExcMessage("Trial did not advance spring history"));
  seam.update_particles(seam.dem_particles, 1., false);
  AssertThrow((seam.get_particle_particle_contact_map()
                 .at(0)
                 .at(1)
                 .tangential_displacement -
               trial_spring)
                  .norm() < 1e-15,
              ExcMessage(
                "Explicit coupling discarded the completed DEM trial"));
  seam.update_particles(seam.dem_particles, 1., true);
  AssertThrow((seam.get_particle_particle_contact_map()
                 .at(0)
                 .at(1)
                 .tangential_displacement -
               accepted_spring)
                  .norm() < 1e-15,
              ExcMessage("Nonlinear trials accumulate contact history"));
  // Separate the pair just beyond contact and check lubrication equivalence.
  seam.dem_particles[0].set_position(0.02, 0);
  seam.dem_particles[1].set_position(0.95, 0);
  interior.dem_particles[1].set_position(-0.05, 0);
  seam.dem_particles[0].velocity[0] = interior.dem_particles[0].velocity[0] =
    -0.1;
  std::vector<Tensor<1, 3>> seam_force(2), interior_force(2), seam_torque(2),
    interior_torque(2);
  seam.calculate_pp_lubrication_force(
    1e-4, 0.02, 0.001, 1., seam_force, seam_torque);
  interior.calculate_pp_lubrication_force(
    1e-4, 0.02, 0.001, 1., interior_force, interior_torque);
  AssertThrow(seam_force[0].norm() > 0. &&
                (seam_force[0] - interior_force[0]).norm() < 1e-12,
              ExcMessage("Periodic lubrication differs"));
  // Two equidistant neighbors must both contribute to lubrication on particle
  // 0.
  auto three_particles = seam.dem_particles;
  auto third           = three_particles[0];
  third.particle_id    = 2;
  third.shape          = third.shape->static_copy();
  third.set_position(0.09, 0);
  third.velocity.clear();
  three_particles.push_back(third);
  IBParticlesDEM<dim> three;
  three.initialize(
    parameters, walls, boundaries, MPI_COMM_WORLD, three_particles);
  three.set_periodic_offset(offset);
  three.update_contact_candidates();
  std::vector<Tensor<1, 3>> three_force(3), three_torque(3);
  three.calculate_pp_lubrication_force(
    1e-4, 0.02, 0.001, 1., three_force, three_torque);
  AssertThrow((three_force[0] - 2. * seam_force[0]).norm() < 1e-12,
              ExcMessage(
                "A lubrication neighbor overwrote another contribution"));
  deallog << dim << "D: normalization, contact, history and lubrication OK"
          << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  initlog();
  test<2>();
  test<3>();
}
