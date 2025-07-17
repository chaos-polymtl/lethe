// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of non-linear (Hertzian) particle-wall
 * contact force is checked.
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
#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>

// Tests (with common definitions)
#include <../tests/dem/test_particles_functions.h>

#include <../tests/tests.h>


using namespace dealii;

template <int dim, typename PropertiesIndex>
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
  const double grid_radius       = 0.5 * GridTools::diameter(tr);
  int          refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  DEMSolverParameters<dim> dem_parameters;
  set_default_dem_parameters(1, dem_parameters);
  auto          &properties = dem_parameters.lagrangian_physical_properties;
  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt                                          = 0.00001;
  double         particle_diameter                           = 0.005;
  unsigned int   rotating_wall_maximum_number                = 6;
  properties.particle_type_number                            = 1;
  properties.youngs_modulus_particle[0]                      = 50000000;
  properties.youngs_modulus_wall                             = 50000000;
  properties.poisson_ratio_particle[0]                       = 0.3;
  properties.poisson_ratio_wall                              = 0.3;
  properties.restitution_coefficient_particle[0]             = 0.5;
  properties.restitution_coefficient_wall                    = 0.5;
  properties.friction_coefficient_particle[0]                = 0.5;
  properties.friction_coefficient_wall                       = 0.5;
  properties.rolling_friction_coefficient_particle[0]        = 0.1;
  properties.rolling_friction_wall                           = 0.1;
  properties.rolling_viscous_damping_coefficient_particle[0] = 0.1;
  properties.rolling_viscous_damping_wall                    = 0.1;
  properties.density_particle[0]                             = 2500;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  // Initializing motion of boundaries
  Tensor<1, dim> translational_and_rotational_velocity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      translational_and_rotational_velocity[d] = 0;
    }
  for (unsigned int counter = 0; counter < rotating_wall_maximum_number;
       ++counter)
    {
      dem_parameters.boundary_conditions.boundary_rotational_speed.insert(
        {counter, 0});
      dem_parameters.boundary_conditions.boundary_translational_velocity.insert(
        {counter, translational_and_rotational_velocity});
      dem_parameters.boundary_conditions.boundary_rotational_vector.insert(
        {counter, translational_and_rotational_velocity});
    }

  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Inserting one particle in contact with a wall
  Point<dim>                       position1 = {-0.998, 0, 0};
  int                              id1       = 0;
  Tensor<1, dim>                   v1{{0.01, 0, 0}};
  Tensor<1, dim>                   omega1{{0, 0, 0}};
  const double                     mass = 1;
  const int                        type = 0;
  Particles::ParticleIterator<dim> pit1 =
    construct_particle_iterator<dim>(particle_handler, tr, position1, id1);

  set_particle_properties<dim, PropertiesIndex>(
    pit1, type, particle_diameter, mass, v1, omega1);

  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  std::vector<double>                          MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  const unsigned int number_of_particles =
    particle_handler.get_max_local_particle_index();
  contact_outcome.resize_interaction_containers(number_of_particles);
  MOI.resize(number_of_particles);
  for (auto &moi_val : MOI)
    moi_val = 1;

  // Finding boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::vector<unsigned int>     outlet_boundaries;
  boundary_cells_object.build(
    tr,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // P-W broad search
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  find_particle_wall_contact_pairs<dim>(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);
  find_particle_wall_contact_pairs<dim>(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // Calling fine search
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_contact_information;
  particle_wall_fine_search<dim>(particle_wall_contact_list,
                                 particle_wall_contact_information);

  // Calling non-linear force
  ParticleWallContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleWallContactForceModel::nonlinear,
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
                           force_object(dem_parameters);
  bool                     floating_wall = false;
  OngoingCollisionLog<dim> ongoing_collision_log;
  CollisionEventLog<dim>   collision_event_log;
  double                   time;
  force_object.calculate_particle_wall_contact(
    particle_wall_contact_information,
    dt,
    time,
    particle_handler,
    floating_wall,
    contact_outcome,
    ongoing_collision_log,
    collision_event_log);

  // Output
  auto particle = particle_handler.begin();
  deallog << "The contact force acting on particle 1 is: "
          << contact_outcome.force[particle->get_id()][0] << " N " << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<3, DEM::DEMProperties::PropertiesIndex>();
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
