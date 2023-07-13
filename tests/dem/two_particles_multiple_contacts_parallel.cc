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

#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
reinitialize_force(Particles::ParticleHandler<dim> &particle_handler,
                   std::vector<Tensor<1, 3>> &      torque,
                   std::vector<Tensor<1, 3>> &      force)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting id of particle as local variable
      unsigned int particle_id = particle->get_id();

      // Reinitializing forces and torques of particles in the system
      force[particle_id][0] = 0;
      force[particle_id][1] = 0;
      force[particle_id][2] = 0;

      torque[particle_id][0] = 0;
      torque[particle_id][1] = 0;
      torque[particle_id][2] = 0;
    }
}

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  double                                    hyper_cube_length = 0.05;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  Tensor<1, 3> g{{0, 0, 0}};
  double       dt                                                    = 0.00001;
  double       particle_diameter                                     = 0.005;
  unsigned int step_end                                              = 1000;
  unsigned int output_frequency                                      = 10;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_particle[0] =
    50000000;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[0] = 0.9;
  dem_parameters.lagrangian_physical_properties
    .restitution_coefficient_particle[0] = 0.9;
  dem_parameters.lagrangian_physical_properties
    .friction_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties
    .rolling_friction_coefficient_particle[0]                       = 0.1;
  dem_parameters.lagrangian_physical_properties.density_particle[0] = 2500;
  double neighborhood_threshold = 1.3 * particle_diameter;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties());

  typename dem_data_structures<2>::particle_index_iterator_map
    local_particle_container;
  typename dem_data_structures<2>::adjacent_particle_pairs
    cleared_local_adjacent_particles;
  typename dem_data_structures<2>::adjacent_particle_pairs
    cleared_ghost_adjacent_particles;

  DEMContactManager<dim> container_manager;

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_object.find_cell_neighbors(
    triangulation,
    container_manager.cells_local_neighbor_list,
    container_manager.cells_ghost_neighbor_list);

  // Creating broad search, fine search and particle-particle force objects
  ParticleParticleBroadSearch<dim> broad_search_object;
  ParticleParticleFineSearch<dim>  fine_search_object;
  ParticleParticleContactForce<
    dim,
    Parameters::Lagrangian::ParticleParticleContactForceModel::
      hertz_mindlin_limit_overlap,
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
                                nonlinear_force_object(dem_parameters);
  VelocityVerletIntegrator<dim> integrator_object;

  MPI_Comm communicator     = triangulation.get_communicator();
  auto     this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

  // Inserting two particles in contact
  Point<2> position1 = {0, 0.007};
  int      id1       = 0;
  Point<2> position2 = {0, 0.001};
  int      id2       = 1;

  // Both particles are inserted the domain owned by process1
  if (this_mpi_process == 1)
    {
      Particles::Particle<dim> particle1(position1, position1, id1);
      typename Triangulation<dim>::active_cell_iterator cell1 =
        GridTools::find_active_cell_around_point(triangulation,
                                                 particle1.get_location());
      Particles::ParticleIterator<dim> pit1 =
        particle_handler.insert_particle(particle1, cell1);
      pit1->get_properties()[DEM::PropertiesIndex::type]    = 0;
      pit1->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
      pit1->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
      pit1->get_properties()[DEM::PropertiesIndex::v_y]     = -0.4;
      pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::mass]    = 1;

      Particles::Particle<dim> particle2(position2, position2, id2);
      typename Triangulation<dim>::active_cell_iterator cell2 =
        GridTools::find_active_cell_around_point(triangulation,
                                                 particle2.get_location());
      Particles::ParticleIterator<dim> pit2 =
        particle_handler.insert_particle(particle2, cell2);
      pit2->get_properties()[DEM::PropertiesIndex::type]    = 0;
      pit2->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
      pit2->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
      pit2->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
      pit2->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::mass]    = 1;
    }


  // Defining variables
  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();

  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (unsigned i = 0; i < MOI.size(); ++i)
    MOI[i] = 1;

  double step_force;

  for (unsigned int iteration = 0; iteration < step_end; ++iteration)
    {
      // Reinitializing forces
      reinitialize_force(particle_handler, torque, force);

      particle_handler.exchange_ghost_particles();

      container_manager.update_local_particles_in_cells(particle_handler,
                                                        false);

      // Calling broad search
      container_manager.execute_particle_particle_broad_search(
        particle_handler);

      // Calling fine search
      container_manager.execute_particle_particle_fine_search(
        neighborhood_threshold);

      // Integration
      // Calling non-linear force
      nonlinear_force_object.calculate_particle_particle_contact_force(
        container_manager, dt, torque, force);

      // Store force before integration for proc 1
      // TODO - Improve this in the future, this is not clean.
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
        step_force = force[0][1];

      // Integration
      integrator_object.integrate(particle_handler, g, dt, torque, force, MOI);

      container_manager.update_contacts();

      if (iteration % output_frequency == 0)
        {
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
            {
              // Output particle 0 position
              for (auto particle = particle_handler.begin();
                   particle != particle_handler.end();
                   ++particle)
                {
                  if (particle->get_id() == 0)
                    {
                      deallog << "The exerted force on particle "
                              << particle->get_id() << " at step " << iteration
                              << " is: " << step_force << std::endl;
                    }
                }
            }
        }
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
