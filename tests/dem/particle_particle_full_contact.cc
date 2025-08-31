// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#include <../tests/dem/full_contact_functions.h>


/**
 * @brief Write time, force, torque and overlap on the log for the full duration of contact.
 * @param output Output of function simulate_full_contact.
 */
void
log_contact_output(full_contact_output &output)
{
  if (output.time.empty())
    {
      deallog << "Error. Contact not found or not ending (too many iterations)."
              << std::endl;
      return;
    }
  const unsigned int time_size = output.time.size();
  for (unsigned int i = 0; i < time_size; ++i)
    {
      deallog << "At time " << std::setw(15) << output.time[i]
              << std::defaultfloat << " force is " << std::setw(10)
              << output.force[i][0] << " " << std::setw(10)
              << output.force[i][1] << " " << std::setw(10)
              << output.force[i][2] << " torque is " << std::setw(10)
              << output.torque[i][0] << " " << std::setw(10)
              << output.torque[i][1] << " " << std::setw(10)
              << output.torque[i][2] << " overlap is " << std::setw(10)
              << output.normal_overlap[i] << std::endl;
    }
}


template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  const int                                 hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  const int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);


  // Defining general simulation parameters
  DEMSolverParameters<dim>                              dem_parameters;
  Parameters::Lagrangian::LagrangianPhysicalProperties &lagrangian_prop =
    dem_parameters.lagrangian_physical_properties;
  Parameters::Lagrangian::ModelParameters<dim> &model_param =
    dem_parameters.model_parameters;

  Tensor<1, dim>     g{{0, 0, 0}};
  const double       dt                                           = 0.00001;
  const double       particle_diameter                            = 0.005;
  const unsigned int output_frequency                             = 10;
  lagrangian_prop.particle_type_number                            = 1;
  lagrangian_prop.youngs_modulus_particle[0]                      = 50000000;
  lagrangian_prop.poisson_ratio_particle[0]                       = 0.3;
  lagrangian_prop.restitution_coefficient_particle[0]             = 0.5;
  lagrangian_prop.friction_coefficient_particle[0]                = 0.5;
  lagrangian_prop.rolling_viscous_damping_coefficient_particle[0] = 0.5;
  lagrangian_prop.rolling_friction_coefficient_particle[0]        = 0.1;
  lagrangian_prop.surface_energy_particle[0]                      = 0.;
  lagrangian_prop.hamaker_constant_particle[0]                    = 0.;
  lagrangian_prop.density_particle[0]                             = 2500;
  model_param.rolling_resistance_method =
    RollingResistanceMethod::constant;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  // Creating containers manager
  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors list
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Setting initial properties of particle 0 and 1
  Point<dim>                       position_0 = {0.4, 0, 0};
  Point<dim>                       position_1 = {0.405, 0, 0};
  Tensor<1, dim>                   velocity_0{{0.01, 0, 0}};
  Tensor<1, dim>                   omega_0{{0, 0, 0}};
  Tensor<1, dim>                   velocity_1{{0, 0, 0}};
  Tensor<1, dim>                   omega_1{{0, 0, 0}};
  initial_particle_properties<dim> particle_properties = {
    {position_0, position_1},              // initial positions
    {0, 1},                                // ids
    {0, 0},                                // types
    {velocity_0, velocity_1},              // initial velocities
    {omega_0, omega_1},                    // initial angular velocities
    {1, 1},                                // masses
    {particle_diameter, particle_diameter} // diameters
  };


  // linear contact force model

  deallog << "For the linear contact force model." << std::endl;

  full_contact_output linear_output =
    simulate_full_contact<dim,
                          PropertiesIndex,
                          ParticleParticleContactForceModel::linear,
                          RollingResistanceMethod::constant>(
      triangulation,
      particle_handler,
      contact_manager,
      dem_parameters,
      particle_properties,
      g,
      dt,
      output_frequency,
      neighborhood_threshold,
      0, // cut_off_factor to allow contact forces without contact
      "linear");

  log_contact_output(linear_output);


  // hertz_mindlin_limit_overlap contact force model

  deallog << " " << std::endl
          << "For the hertz_mindlin_limit_overlap contact force model."
          << std::endl;

  full_contact_output hmlo_output = simulate_full_contact<
    dim,
    PropertiesIndex,
    ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
    RollingResistanceMethod::constant>(triangulation,
                                                  particle_handler,
                                                  contact_manager,
                                                  dem_parameters,
                                                  particle_properties,
                                                  g,
                                                  dt,
                                                  output_frequency,
                                                  neighborhood_threshold,
                                                  0,
                                                  "hmlo");

  log_contact_output(hmlo_output);
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
