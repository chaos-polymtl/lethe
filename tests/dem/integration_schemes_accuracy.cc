// SPDX-FileCopyrightText: Copyright (c) 2021-2022, 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, accuracies of different integration
 * schemes are compared using a pendulum oscillation model
 */

// Deal.II includes
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

#include <dem/explicit_euler_integrator.h>
#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)
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
  int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining simulation general parameters
  Tensor<1, dim>     g{{0, 0, 0}};
  double             dt1             = 0.1;
  double             dt2             = 0.05;
  const unsigned int time_step_ratio = dt1 / dt2;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  int    id              = 0;
  double particle_mass   = 1;
  double spring_constant = 1;


  DEMSolverParameters<dim> dem_parameters;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;

  // Initial condition
  double   t         = 0;
  double   x0        = 0.3;
  Point<3> position1 = {0, 0, x0};
  double   t_final   = 0.999999;
  double   x_analytical;
  double   particle_axial_position_error_Euler_dt1;
  double   particle_axial_position_error_Euler_dt2;
  double   particle_axial_position_error_Verlet_dt1;
  double   particle_axial_position_error_Verlet_dt2;

  Particles::Particle<dim> particle0(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle0_cell =
    GridTools::find_active_cell_around_point(tr, particle0.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit0 =
    particle_handler.insert_particle(particle0, particle0_cell);

  pit0->get_properties()[DEM::DEMProperties::PropertiesIndex::v_x] = 0;
  pit0->get_properties()[DEM::DEMProperties::PropertiesIndex::v_y] = 0;
  pit0->get_properties()[DEM::DEMProperties::PropertiesIndex::v_z] = 0;
  pit0->get_properties()[DEM::DEMProperties::PropertiesIndex::mass] =
    particle_mass;

  // Calling integrators
  ExplicitEulerIntegrator<dim, DEM::DEMProperties::PropertiesIndex>
    explicit_euler_object;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  MOI[0] = 1.;

  // Explicit Euler
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      while (t < t_final)
        {
          Tensor<1, dim> force_tensor;
          force_tensor[dim - 1] =
            -spring_constant * particle_iterator->get_location()[dim - 1];
          force[particle_iterator->get_local_index()] = force_tensor;
          explicit_euler_object.integrate(
            particle_handler, g, dt1, torque, force, MOI);

          t += dt1;
        }
      // Output Analytical
      x_analytical = x0 * cos(sqrt(spring_constant / particle_mass) * (t));
      particle_axial_position_error_Euler_dt1 =
        particle_iterator->get_location()[dim - 1] - x_analytical;
    }

  particle_handler.clear_particles();
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle1_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle1_cell);


  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_x] = 0;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_y] = 0;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_z] = 0;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::mass] =
    particle_mass;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());

  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      t = 0;

      while (t < t_final)
        {
          Tensor<1, dim> force_tensor;
          force_tensor[dim - 1] =
            -spring_constant * particle_iterator->get_location()[dim - 1];
          force[particle_iterator->get_local_index()] = force_tensor;
          explicit_euler_object.integrate(
            particle_handler, g, dt2, torque, force, MOI);
          t += dt2;
        }
      // Output Analytical
      x_analytical = x0 * cos(sqrt(spring_constant / particle_mass) * (t));
      particle_axial_position_error_Euler_dt2 =
        particle_iterator->get_location()[dim - 1] - x_analytical;
    }
  deallog << "Explicit Euler is a "
          << std::abs(particle_axial_position_error_Euler_dt1 /
                      particle_axial_position_error_Euler_dt2) /
               time_step_ratio
          << " order integration scheme" << std::endl;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle2(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle2_cell =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, particle2_cell);

  pit2->get_properties()[DEM::DEMProperties::PropertiesIndex::v_x] = 0;
  pit2->get_properties()[DEM::DEMProperties::PropertiesIndex::v_y] = 0;
  pit2->get_properties()[DEM::DEMProperties::PropertiesIndex::v_z] = 0;
  pit2->get_properties()[DEM::DEMProperties::PropertiesIndex::mass] =
    particle_mass;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());

  // Create Velocity Verlet integrator
  VelocityVerletIntegrator<dim, DEM::DEMProperties::PropertiesIndex>
    velocity_verlet_object;

  // Output Velocity Verlet
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      t = 0;

      force[particle_iterator->get_local_index()][dim - 1] = -x0;

      velocity_verlet_object.integrate_half_step_location(
        particle_handler, g, dt1, torque, force, MOI);
      t += dt1;

      while (t < t_final)
        {
          Tensor<1, dim> force_tensor;
          force_tensor[dim - 1] =
            -spring_constant * particle_iterator->get_location()[dim - 1];
          force[particle_iterator->get_local_index()] = force_tensor;
          velocity_verlet_object.integrate(
            particle_handler, g, dt1, torque, force, MOI);

          t += dt1;
        }
      // Output Analytical
      x_analytical = x0 * cos(sqrt(spring_constant / particle_mass) * (t));
      particle_axial_position_error_Verlet_dt1 =
        particle_iterator->get_location()[dim - 1] - x_analytical;
    }
  particle_handler.clear_particles();
  Particles::Particle<dim> particle3(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle3_cell =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, particle3_cell);

  pit3->get_properties()[DEM::DEMProperties::PropertiesIndex::v_x] = 0;
  pit3->get_properties()[DEM::DEMProperties::PropertiesIndex::v_y] = 0;
  pit3->get_properties()[DEM::DEMProperties::PropertiesIndex::v_z] = 0;
  pit3->get_properties()[DEM::DEMProperties::PropertiesIndex::mass] =
    particle_mass;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());

  // Output Velocity Verlet
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      t = 0;


      force[particle_iterator->get_local_index()][dim - 1] = -x0;
      velocity_verlet_object.integrate_half_step_location(
        particle_handler, g, dt2, torque, force, MOI);
      t += dt2;

      while (t < t_final)
        {
          Tensor<1, dim> force_tensor;
          force_tensor[dim - 1] =
            -spring_constant * particle_iterator->get_location()[dim - 1];
          force[particle_iterator->get_local_index()] = force_tensor;

          velocity_verlet_object.integrate(
            particle_handler, g, dt2, torque, force, MOI);
          t += dt2;
        }
      // Output Analytical
      x_analytical = x0 * cos(sqrt(spring_constant / particle_mass) * (t));
      particle_axial_position_error_Verlet_dt2 =
        particle_iterator->get_location()[dim - 1] - x_analytical;
    }

  deallog << "Velocity Verlet is a "
          << std::abs(particle_axial_position_error_Verlet_dt1 /
                      particle_axial_position_error_Verlet_dt2) /
               time_step_ratio
          << " order integration scheme" << std::endl;
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
