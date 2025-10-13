// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test if the particle projector gives a valid solution for a constant void fraction.
 * The constant void fraction is defined by homogenously placing particles in a
 * grid within a cubic triangulation and averaging at a significantly larger
 * volume than the volume of the particles. Right now, the void fraction is not
 * that homogenous due to the effect of the walls on the QCM averaging volume
 */

// Deal.II includes
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
// Lethe
#include <core/dem_properties.h>

#include <fem-dem/particle_projector.h>

// Tests
#include <../tests/tests.h>
#include <../tests/tests_utilities.h>

using namespace dealii;


// Fill an empty particle handler that is ready for CFD-DEM simulation using a
// list of points

template <int dim>
void
generate_particle_grid(const Point<dim>          pt1,
                       const Point<dim>          pt2,
                       const unsigned int        n_refinements,
                       const Triangulation<dim> &background_triangulation,
                       Particles::ParticleHandler<dim> &particle_handler)
{
  parallel::distributed::Triangulation<dim> particle_triangulation(
    MPI_COMM_WORLD);


  GridGenerator::hyper_rectangle(particle_triangulation, pt1, pt2, false);
  particle_triangulation.refine_global(n_refinements);

  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    background_triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

  std::vector<std::vector<double>> properties(
    particle_triangulation.n_locally_owned_active_cells(),
    std::vector<double>(DEM::CFDDEMProperties::n_properties, 0.));

  const MappingQ1<dim> mapping;

  Particles::Generators::quadrature_points(particle_triangulation,
                                           QMidpoint<dim>(),
                                           global_bounding_boxes,
                                           particle_handler,
                                           mapping,
                                           properties);

  deallog << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
}

void
test_void_fraction_qcm(const unsigned int fe_degree,
                       const unsigned int number_quadrature_points)
{
  const auto         my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  constexpr unsigned int n_refinements = 1;

  // We make a background triangulation which consists in a 1x1x1 cube.
  parallel::distributed::Triangulation<3> domain_triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(domain_triangulation, 0, 1, false);

  // We refine this triangulation for a given number of refinement
  domain_triangulation.refine_global(n_refinements);

  // We initialze the particle handler
  const MappingQ1<3>            mapping;
  Particles::ParticleHandler<3> particle_handler;
  particle_handler.initialize(domain_triangulation,
                              mapping,
                              DEM::CFDDEMProperties::n_properties);


  // We generate the grid of particles
  generate_particle_grid(Point<3>({0, 0, 0}),
                         Point<3>({1, 1, 1}),
                         5,
                         domain_triangulation,
                         particle_handler);

  // We fix the diameter of all the particles to have a void fraction of 0.9 and
  // a solid fraction of 0.1. The volume of the underlying grid is 1.
  const double target_eps_solid = 0.1;
  const double dp =
    std::pow(target_eps_solid * 1. / particle_handler.n_global_particles() * 6 /
               numbers::PI,
             1. / 3.);

  deallog << "dp = " << dp << std::endl;

  // Loop over all the particles and set their diameter to dp
  for (auto &particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      particle_properties[DEM::CFDDEMProperties::dp] = dp;
    }

  // Setup the ParticleProjector. For this we need VoidFractionParameters,
  // LinearSolverParameters and a conditional OSStream.

  // Setup default VoidFractionParameters
  std::shared_ptr<Parameters::VoidFractionParameters<3>>
    void_fraction_parameters = make_default_void_fraction_parameters();
  void_fraction_parameters->n_quadrature_points = number_quadrature_points;

  // Setup a default linear solver.
  Parameters::LinearSolver linear_solver_parameters =
    make_default_linear_solver();

  // Setup a pcout which is required by the ParticleProjector.

  ConditionalOStream pcout(deallog.get_file_stream(),
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  BoundaryConditions::NSBoundaryConditions<3> boundary_conditions;

  ParticleProjector<3> particle_projector(&domain_triangulation,
                                          void_fraction_parameters,
                                          linear_solver_parameters,
                                          &particle_handler,
                                          fe_degree,
                                          false,
                                          pcout);



  // The particle projector requires that its DOF values be set-up and the
  // constraints as well.
  particle_projector.setup_dofs();
  particle_projector.setup_constraints(boundary_conditions);

  // Calculate the void fraction and print the DOF values
  particle_projector.calculate_void_fraction(0.);

  for (unsigned int r = 0; r < n_procs; ++r)
    {
      if (my_rank == r)
        {
          deallog << "Rank " << r << " owns: ";
          for (const auto &i : particle_projector.void_fraction_solution)
            deallog << i << " ";
          deallog << std::endl;
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }


  // Integrate the void fraction field over the cells to ensure void fraction
  // conservation
  FEValues<3> fe_values_void_fraction(*particle_projector.mapping,
                                      *particle_projector.fe,
                                      *particle_projector.quadrature,
                                      update_values | update_JxW_values);

  double              total_particle_volume = 0;
  std::vector<double> void_fraction_values(
    fe_values_void_fraction.n_quadrature_points);
  for (const auto &cell :
       particle_projector.dof_handler.active_cell_iterators())
    {
      fe_values_void_fraction.reinit((cell));
      fe_values_void_fraction.get_function_values(
        particle_projector.void_fraction_solution, void_fraction_values);
      for (unsigned int q = 0; q < fe_values_void_fraction.n_quadrature_points;
           ++q)
        total_particle_volume += fe_values_void_fraction.get_JxW_values()[q] *
                                 (1. - void_fraction_values[q]);
    }
  deallog << "Total particle volume " << total_particle_volume << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      deallog << "Void fraction: fe_degree=1    number_quadrature_points=2"
              << std::endl;
      test_void_fraction_qcm(1, 2);
      deallog << "Void fraction: fe_degree=1    number_quadrature_points=3"
              << std::endl;
      test_void_fraction_qcm(1, 3);
      deallog << "Void fraction: fe_degree=2    number_quadrature_points=3"
              << std::endl;
      test_void_fraction_qcm(2, 3);
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
