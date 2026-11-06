// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the QCM projection of the particle temperature on a mesh with
 * hanging nodes.
 *
 * The particles carry the temperature T = 300 + 10 x. The test prints the
 * integral of the projected temperature provided to the multiphysics interface
 * and its L2 distance to T. These are regression values: the projection is not
 * exactly T, since the QCM filter is truncated near the boundary and the
 * smoothing imposes a zero normal gradient there.
 */

// Deal.II includes
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/boundary_conditions.h>
#include <core/dem_properties.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>

#include <fem-dem/particle_projector.h>

#include <cmath>
#include <exception>
#include <iomanip>
#include <memory>
#include <vector>

// Tests
#include <../tests/tests.h>
#include <../tests/tests_utilities.h>

using namespace dealii;

/**
 * @brief Fill an empty particle handler with one particle at the center of
 * each cell of a uniformly refined grid of the box [pt1, pt2].
 *
 * @param[in] pt1 Lower corner of the box of particles.
 * @param[in] pt2 Upper corner of the box of particles.
 * @param[in] n_refinements Number of global refinements of the particle grid.
 * @param[in] background_triangulation Triangulation of the fluid domain.
 * @param[in,out] particle_handler Particle handler that receives the particles.
 */
void
generate_particle_grid(const Point<3>                &pt1,
                       const Point<3>                &pt2,
                       const unsigned int             n_refinements,
                       const Triangulation<3>        &background_triangulation,
                       Particles::ParticleHandler<3> &particle_handler)
{
  parallel::distributed::Triangulation<3> particle_triangulation(
    MPI_COMM_WORLD);

  GridGenerator::hyper_rectangle(particle_triangulation, pt1, pt2, false);
  particle_triangulation.refine_global(n_refinements);

  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    background_triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  std::vector<std::vector<double>> properties(
    particle_triangulation.n_locally_owned_active_cells(),
    std::vector<double>(DEM::CFDDEMMPProperties::n_properties, 0.));

  const MappingQ1<3> mapping;

  Particles::Generators::quadrature_points(particle_triangulation,
                                           QMidpoint<3>(),
                                           global_bounding_boxes,
                                           particle_handler,
                                           mapping,
                                           properties);

  deallog << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
}

/**
 * @brief Project the particle temperature T = 300 + 10 x, then print the
 * integral of the projection and its L2 distance to T.
 */
void
test()
{
  // Particle temperature T = temperature_at_origin + temperature_gradient * x
  const double temperature_at_origin = 300.;
  const double temperature_gradient  = 10.;

  // Unit cube refined twice, then once more in the octant [0, 0.5]^3, which
  // creates hanging nodes at the interface between the two refinement levels.
  parallel::distributed::Triangulation<3> domain_triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(domain_triangulation, 0., 1., false);
  domain_triangulation.refine_global(2);

  for (const auto &cell : domain_triangulation.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      const Point<3> center = cell->center();
      if (center[0] < 0.5 && center[1] < 0.5 && center[2] < 0.5)
        cell->set_refine_flag();
    }
  domain_triangulation.execute_coarsening_and_refinement();

  deallog << "Number of cells: " << domain_triangulation.n_global_active_cells()
          << std::endl;

  // Particles with thermal properties, one at the center of each cell of a
  // 16 x 16 x 16 grid
  const MappingQ1<3>            mapping;
  Particles::ParticleHandler<3> particle_handler;
  particle_handler.initialize(domain_triangulation,
                              mapping,
                              DEM::CFDDEMMPProperties::n_properties);

  generate_particle_grid(Point<3>({0., 0., 0.}),
                         Point<3>({1., 1., 1.}),
                         4,
                         domain_triangulation,
                         particle_handler);

  // Particle diameter that gives a solid fraction of 0.1 in the unit cube
  const double target_solid_fraction = 0.1;
  const double particle_diameter =
    std::pow(target_solid_fraction / particle_handler.n_global_particles() *
               6. / numbers::PI,
             1. / 3.);

  for (auto particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      particle_properties[DEM::CFDDEMMPProperties::dp] = particle_diameter;
      particle_properties[DEM::CFDDEMMPProperties::T] =
        temperature_at_origin +
        temperature_gradient * particle.get_location()[0];
    }

  // Exchange the ghost particles once their properties are set
  particle_handler.exchange_ghost_particles(true);

  // The projection of the particle temperature requires an l2 smoothing
  // length > 0.
  std::shared_ptr<Parameters::VoidFractionParameters<3>>
    void_fraction_parameters     = make_default_void_fraction_parameters();
  void_fraction_parameters->mode = Parameters::VoidFractionMode::qcm;
  void_fraction_parameters->l2_smoothing_length = 0.1;

  Parameters::LinearSolver linear_solver_parameters =
    make_default_linear_solver();
  linear_solver_parameters.minimum_residual = 1e-12;
  linear_solver_parameters.verbosity        = Parameters::Verbosity::quiet;

  ConditionalOStream pcout(deallog.get_file_stream(), true);

  BoundaryConditions::NSBoundaryConditions<3> boundary_conditions;

  ParticleProjector<3, DEM::CFDDEMMPProperties::PropertiesIndex>
    particle_projector(&domain_triangulation,
                       void_fraction_parameters,
                       linear_solver_parameters,
                       &particle_handler,
                       1,
                       false,
                       pcout);

  particle_projector.setup_dofs();
  particle_projector.setup_constraints(boundary_conditions);

  particle_projector.calculate_void_fraction(0.);
  particle_projector.calculate_particle_temperature_projection();

  // Integral of the projected temperature provided to the multiphysics
  // interface, and its L2 distance to T
  const auto     &temperature_field = particle_projector.particle_temperature;
  const QGauss<3> quadrature(temperature_field.fe->degree + 1);
  FEValues<3>     fe_values(*particle_projector.mapping,
                        *temperature_field.fe,
                        quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values);

  std::vector<double> projected_temperature(quadrature.size());
  double              integral            = 0.;
  double              squared_l2_distance = 0.;

  for (const auto &cell :
       temperature_field.dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_values(
        *particle_projector.particle_temperature_locally_relevant,
        projected_temperature);

      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          const double particle_temperature =
            temperature_at_origin +
            temperature_gradient * fe_values.quadrature_point(q)[0];

          integral += projected_temperature[q] * fe_values.JxW(q);
          squared_l2_distance +=
            Utilities::fixed_power<2>(projected_temperature[q] -
                                      particle_temperature) *
            fe_values.JxW(q);
        }
    }

  deallog << "Integral of the projected temperature: "
          << Utilities::MPI::sum(integral, MPI_COMM_WORLD) << std::endl;
  deallog << "L2 distance to 300 + 10 x: "
          << std::sqrt(Utilities::MPI::sum(squared_l2_distance, MPI_COMM_WORLD))
          << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;
      deallog << std::fixed << std::setprecision(3);

      test();
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
