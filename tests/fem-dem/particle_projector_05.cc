// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the QCM projection of a field (here the particle temperature is
 * used for the test case) on a mesh with hanging nodes.
 *
 * All the particles carry the same temperature T_0. The QCM average of the
 * particle temperature is therefore T_0 at every quadrature point that sees
 * particles, and the constant field T_0 solves the smoothed L2 projection
 * exactly, including at the hanging nodes, since the constrained finite
 * element space contains the constant functions. The test prints the minimum
 * and the maximum of the projected temperature over the degrees of freedom,
 * which must both be T_0:
 * - in the deal.II vector of the projection, on the locally owned cells;
 * - in its copy provided to the multiphysics interface, on the locally owned
 *   and ghost cells, which also checks the ghost values of the copy.
 */

// Deal.II includes
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/vector.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/boundary_conditions.h>
#include <core/dem_properties.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>

#include <fem-dem/particle_projector.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
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

void
test()
{
  // Temperature carried by all the particles
  const double particle_temperature = 300.;

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
  // 16 x 16 x 16 grid (the grid for the particles is independent of the grid
  // for the fluid domain and is uniformly refined everywhere).
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

  // Exchange ghost particles after particle_triangulation (inside
  // generate_particle_grid) has gone out of scope so that its destructor does
  // not overlap with MPI operations from the ghost exchange.
  particle_handler.exchange_ghost_particles(true);

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
      particle_properties[DEM::CFDDEMMPProperties::T]  = particle_temperature;
    }

  // The projection of the particle temperature requires the qcm mode and an
  // l2 smoothing length > 0.
  std::shared_ptr<Parameters::VoidFractionParameters<3>>
    void_fraction_parameters = make_default_void_fraction_parameters();
  void_fraction_parameters->l2_smoothing_length          = 0.1;
  void_fraction_parameters->project_particle_temperature = true;

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

  // The particle temperature is projected, and copied for the multiphysics
  // interface, right after the void fraction.
  particle_projector.calculate_void_fraction(0.);

  const auto &temperature_field = particle_projector.particle_temperature;

  // Number of degrees of freedom constrained by the hanging nodes
  unsigned int n_constrained_dofs = 0;
  for (const auto dof : temperature_field.locally_owned_dofs)
    {
      if (temperature_field.particle_field_constraints.is_constrained(dof))
        ++n_constrained_dofs;
    }
  n_constrained_dofs = Utilities::MPI::sum(n_constrained_dofs, MPI_COMM_WORLD);

  deallog << "Number of constrained temperature DoFs: " << n_constrained_dofs
          << std::endl;

  // Minimum and maximum of the projected temperature, for the deal.II vector
  // of the projection (locally owned cells) and for its copy in the interface
  // (locally owned and ghost cells).
  Vector<double> cell_values(temperature_field.fe->n_dofs_per_cell());
  double         min_projection     = std::numeric_limits<double>::max();
  double         max_projection     = std::numeric_limits<double>::lowest();
  double         min_interface_copy = std::numeric_limits<double>::max();
  double         max_interface_copy = std::numeric_limits<double>::lowest();

  for (const auto &cell :
       temperature_field.dof_handler->active_cell_iterators())
    {
      if (cell->is_artificial())
        continue;

      if (cell->is_locally_owned())
        {
          cell->get_dof_values(temperature_field.particle_field_solution,
                               cell_values);
          for (const double value : cell_values)
            {
              min_projection = std::min(min_projection, value);
              max_projection = std::max(max_projection, value);
            }
        }

      cell->get_dof_values(
        *temperature_field.particle_field_locally_relevant_solution,
        cell_values);
      for (const double value : cell_values)
        {
          min_interface_copy = std::min(min_interface_copy, value);
          max_interface_copy = std::max(max_interface_copy, value);
        }
    }

  min_projection     = Utilities::MPI::min(min_projection, MPI_COMM_WORLD);
  max_projection     = Utilities::MPI::max(max_projection, MPI_COMM_WORLD);
  min_interface_copy = Utilities::MPI::min(min_interface_copy, MPI_COMM_WORLD);
  max_interface_copy = Utilities::MPI::max(max_interface_copy, MPI_COMM_WORLD);

  deallog << "Projected temperature: min = " << min_projection
          << ", max = " << max_projection << " (expected "
          << particle_temperature << ")" << std::endl;
  deallog << "Interface copy, ghosts included: min = " << min_interface_copy
          << ", max = " << max_interface_copy << " (expected "
          << particle_temperature << ")" << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;

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
