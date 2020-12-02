/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

/**
 * @brief This test reports the normal overlap and corresponding normal force
 * during a complete particle-wall contact. Interested reader may be interested
 * inplotting normal force against normal overlap.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_contact_force.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  std::vector<unsigned int> holes;
  holes = {1, 1};
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::cheese(tr, holes);
  int refinement_number = 1;
  tr.refine_global(refinement_number);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;
  unsigned int             step              = 0;
  const int                writing_frequency = 100;

  // Defining general simulation parameters
  const unsigned int n_properties = 21;
  Tensor<1, dim>     g{{0, -9.81}};
  double             dt                                         = 0.0001;
  double             particle_diameter                          = 0.1;
  int                particle_density                           = 2000;
  dem_parameters.physical_properties.particle_type_number       = 1;
  dem_parameters.physical_properties.youngs_modulus_particle[0] = 200000000000;
  dem_parameters.physical_properties.youngs_modulus_wall        = 200000000000;
  dem_parameters.physical_properties.poisson_ratio_particle[0]  = 0.3;
  dem_parameters.physical_properties.poisson_ratio_wall         = 0.3;
  dem_parameters.physical_properties.restitution_coefficient_particle[0] = 0.95;
  dem_parameters.physical_properties.restitution_coefficient_wall        = 0.95;
  dem_parameters.physical_properties.friction_coefficient_particle[0]    = 0.05;
  dem_parameters.physical_properties.friction_coefficient_wall           = 0.05;
  dem_parameters.physical_properties.rolling_friction_coefficient_particle[0] =
    0.1;
  dem_parameters.physical_properties.rolling_friction_wall = 0.1;
  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  // Inserting one particle in contact with wall
  Point<dim>               position1 = {0.97, 2.05};
  int                      id        = 0;
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle_cell);
  pit1->get_properties()[DEM::PropertiesIndex::type]        = 0;
  pit1->get_properties()[DEM::PropertiesIndex::dp]          = particle_diameter;
  pit1->get_properties()[DEM::PropertiesIndex::rho]         = particle_density;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass]        = 1;
  pit1->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Construct boundary cells object and build it
  BoundaryCellsInformation<dim> boundary_cells_object;
  boundary_cells_object.build(tr);

  // Particle-point broad search
  std::unordered_map<int,
                     std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
                                    contact_candidates;
  ParticlePointLineBroadSearch<dim> broad_search_object;

  // Particle-point fine search
  ParticlePointLineFineSearch<dim> fine_search_object;
  std::unordered_map<int, particle_point_line_contact_info_struct<dim>>
    contact_information;

  ParticlePointLineForce<dim>   force_object;
  VelocityVerletIntegrator<dim> integrator_object;

  for (double time = 0; time < 0.2; time += dt)
    {
      auto particle = particle_handler.begin();
      particle->get_properties()[DEM::PropertiesIndex::force_x] = 0;
      particle->get_properties()[DEM::PropertiesIndex::force_y] = 0;

      contact_candidates =
        broad_search_object.find_particle_point_contact_pairs(
          particle_handler,
          boundary_cells_object.get_boundary_cells_with_points());

      contact_information =
        fine_search_object.particle_point_fine_search(contact_candidates,
                                                      neighborhood_threshold);

      force_object.calculate_particle_point_contact_force(
        &contact_information, dem_parameters.physical_properties);
      integrator_object.integrate(particle_handler, g, dt);

      if (step % writing_frequency == 0)
        {
          deallog << particle->get_location() << std::endl;
        }
      ++step;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<2>();
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
