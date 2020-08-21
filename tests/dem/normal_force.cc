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

// This test reports the normal overlap and corresponding normal force during a
// complete particle-wall contact. Interested reader may be interested in
// plotting normal force against normal overlap

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_force.h>
#include <dem/pw_fine_search.h>
#include <dem/pw_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  const unsigned int n_properties = 21;
  Tensor<1, dim>     g{{0, 0, 0}};
  double             dt                                      = 0.000001;
  double             particle_diameter                       = 0.001;
  int                particle_density                        = 7850;
  dem_parameters.physical_properties.Youngs_modulus_particle = 200000000000;
  dem_parameters.physical_properties.Youngs_modulus_wall     = 200000000000;
  dem_parameters.physical_properties.Poisson_ratio_particle  = 0.3;
  dem_parameters.physical_properties.Poisson_ratio_wall      = 0.3;
  dem_parameters.physical_properties.restitution_coefficient_particle = 0.5;
  dem_parameters.physical_properties.restitution_coefficient_wall     = 0.5;
  dem_parameters.physical_properties.friction_coefficient_particle    = 0.3;
  dem_parameters.physical_properties.friction_coefficient_wall        = 0.3;
  dem_parameters.physical_properties.rolling_friction_particle        = 0.1;
  dem_parameters.physical_properties.rolling_friction_wall            = 0.1;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  // Inserting one particle in contact with wall
  Point<dim>               position1 = {-0.999, 0, 0};
  int                      id        = 0;
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle_cell);
  pit1->get_properties()[0]  = id;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = particle_diameter;
  pit1->get_properties()[3]  = particle_density;
  pit1->get_properties()[4]  = -1.0;
  pit1->get_properties()[5]  = 0;
  pit1->get_properties()[6]  = 0;
  pit1->get_properties()[7]  = 0;
  pit1->get_properties()[8]  = 0;
  pit1->get_properties()[9]  = 0;
  pit1->get_properties()[10] = 0;
  pit1->get_properties()[11] = 0;
  pit1->get_properties()[12] = 0;
  pit1->get_properties()[13] = 0;
  pit1->get_properties()[14] = 0;
  pit1->get_properties()[15] = 0;
  pit1->get_properties()[16] = 1;
  pit1->get_properties()[17] = 1;

  // Finding boundary cells
  std::vector<typename Triangulation<dim>::active_cell_iterator>
                                                 boundary_cells_with_faces;
  std::map<int, boundary_cells_info_struct<dim>> boundary_cell_information;
  FindBoundaryCellsInformation<dim>              boundary_cells_object;
  boundary_cell_information =
    boundary_cells_object.find_boundary_cells_information(
      boundary_cells_with_faces, tr);

  // P-W broad search
  PWBroadSearch<dim> pw_broad_search_object;
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    pw_contact_list;
  pw_broad_search_object.find_PW_Contact_Pairs(boundary_cell_information,
                                               particle_handler,
                                               pw_contact_list);

  // P-W fine search
  PWFineSearch<dim> pw_fine_search_object;
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
                                pw_contact_information;
  PWNonLinearForce<dim>         pw_force_object;
  VelocityVerletIntegrator<dim> integrator_object;
  double                        distance;

  for (double time = 0; time < 0.00115; time += dt)
    {
      auto particle = particle_handler.begin();
      particle->get_properties()[DEM::PropertiesIndex::force_x] = 0;
      particle->get_properties()[DEM::PropertiesIndex::force_y] = 0;
      if (dim == 3)
        {
          particle->get_properties()[DEM::PropertiesIndex::force_z] = 0;
        }
      distance = hyper_cube_length + particle->get_location()[0] -
                 particle->get_properties()[DEM::PropertiesIndex::dp] / 2.0;

      if (distance > 0.0)
        {
          // If particle and wall are not in contact, only the integration class
          // is called
          integrator_object.integrate(particle_handler, g, dt);
        }
      else
        {
          // If particle and wall are in contact
          pw_fine_search_object.pw_Fine_Search(pw_contact_list,
                                               pw_contact_information);
          auto pw_pairs_in_contact_iterator =
            &pw_contact_information.begin()->second;
          auto pw_contact_information_iterator =
            pw_pairs_in_contact_iterator->begin();

          pw_contact_information_iterator->second.tangential_overlap[0] = 0.0;
          pw_contact_information_iterator->second.tangential_overlap[1] = 0.0;
          if (dim == 3)
            {
              pw_contact_information_iterator->second.tangential_overlap[2] =
                0.0;
            }
          pw_contact_information_iterator->second
            .tangential_relative_velocity[0] = 0.0;
          pw_contact_information_iterator->second
            .tangential_relative_velocity[1] = 0.0;
          if (dim == 3)
            {
              pw_contact_information_iterator->second
                .tangential_relative_velocity[2] = 0.0;
            }

          pw_force_object.calculate_pw_contact_force(&pw_contact_information,
                                                     dem_parameters,
                                                     dt);
          integrator_object.integrate(particle_handler, g, dt);

          deallog << " "
                  << pw_contact_information_iterator->second.normal_overlap
                  << " "
                  << particle->get_properties()[DEM::PropertiesIndex::force_x]
                  << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
