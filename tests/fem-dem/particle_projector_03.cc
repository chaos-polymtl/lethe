// SPDX-FileCopyrightText: Copyright (c) 2015-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test if the particle projector gives a valid solution for a constant
 * and a linear forcing field stored onto the particles.
 *
 * This test is very similar to particle_projector_02 with two main distinction:
 * 1. The mesh for the domain is refined further to ensure that every
 *    processor owns some cells.
 * 2. The test is adapted to be run in parallel. Essentially, it tests
 * furthermore the projection capabilities in both serial and parallel.
 *
 * The test consists of the following steps (which are the same as
 * projector_02). First, the test is run with a constant forcing applied to the
 * particles. The goal of the test is to verify if the forcing is conservative.
 * Afterward, the test is run with a linear force function. The goal of the
 * test is then to ensure that the forcing is conservative, but also that
 * it preserves the function that originally described the forcing (in this
 * case a linear profile).
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

  particle_handler.exchange_ghost_particles(true);

  deallog << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
}

template <int dim>
class LinearForce : public Function<dim>
{
public:
  LinearForce()
    : Function<dim>(dim)
  {}
  virtual double
  value(const Point<dim> &p, const unsigned int component) const override
  {
    if (component == 0)
      return p[1];
    if (component == 1)
      return p[2];
    else
      return p[0];
  }
};

void
test_void_fraction_qcm(const unsigned int fe_degree,
                       const unsigned int number_quadrature_points,
                       Function<3>       &force_distribution,
                       const bool         output_vtu,
                       const std::string &vtu_label)
{
  const auto         my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  constexpr unsigned int n_refinements = 2;

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
    std::pow(target_eps_solid * 1. / particle_handler.n_global_particles() * 3 *
               8. / 4. / numbers::PI,
             1. / 3.);

  deallog << "dp = " << dp << std::endl;

  Tensor<1, 3> total_particle_force_on_particles({0., 0., 0.});
  // Loop over all the particles and set their diameter to dp
  // and the force to a constant value
  for (auto particle : particle_handler)
    {
      auto particle_properties = particle.get_properties();
      particle_properties[DEM::CFDDEMProperties::dp] = dp;
      for (unsigned int d = 0; d < 3; ++d)
        {
          particle_properties
            [DEM::CFDDEMProperties::fem_force_two_way_coupling_x + d] =
              force_distribution.value(particle.get_location(), d);
          total_particle_force_on_particles[d] += particle_properties
            [DEM::CFDDEMProperties::fem_force_two_way_coupling_x + d];
        }
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

  linear_solver_parameters.minimum_residual = 1e-10;
  linear_solver_parameters.verbosity        = Parameters::Verbosity::quiet;

  // Setup a pcout which is required by the ParticleProjector.

  ConditionalOStream pcout(deallog.get_file_stream(), true);

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

  // Calculate the void fraction to initialize the QCM weights
  particle_projector.calculate_void_fraction(0.);

  // Calculate the force projection
  // For this we only need to calculate the field projection!
  announce_string(pcout, "Force on particle");
  particle_projector.calculate_field_projection(
    particle_projector.fluid_force_on_particles_two_way_coupling);


  for (unsigned int r = 0; r < n_procs; ++r)
    {
      if (my_rank == r)
        {
          deallog << "Rank " << r << " owns: ";
          for (const auto &i :
               particle_projector.fluid_force_on_particles_two_way_coupling
                 .particle_field_solution)
            deallog << i << " ";
          deallog << std::endl;
        }
    }

  // Integrate the force field over the cells to check force conservation
  FEValues<3> fe_values(
    *particle_projector.mapping,
    *particle_projector.fluid_force_on_particles_two_way_coupling.fe,
    *particle_projector.quadrature,
    update_values | update_JxW_values);

  // Field extractor if dim components are used
  FEValuesExtractors::Vector vector_extractor;
  vector_extractor.first_vector_component = 0;
  Tensor<1, 3>              total_particle_force_on_fluid({0, 0, 0});
  std::vector<Tensor<1, 3>> force_values(fe_values.n_quadrature_points);
  for (const auto &cell :
       particle_projector.fluid_force_on_particles_two_way_coupling.dof_handler
         .active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit((cell));
          fe_values[vector_extractor].get_function_values(
            particle_projector.fluid_force_on_particles_two_way_coupling
              .particle_field_solution,
            force_values);
          for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
            total_particle_force_on_fluid +=
              fe_values.get_JxW_values()[q] * force_values[q];
        }
    }

  // Sum the total force across all the processors
  total_particle_force_on_fluid =
    Utilities::MPI::sum(total_particle_force_on_fluid,
                        domain_triangulation.get_mpi_communicator());
  total_particle_force_on_particles =
    Utilities::MPI::sum(total_particle_force_on_particles,
                        domain_triangulation.get_mpi_communicator());

  deallog << "Total particle force on the fluid mesh "
          << total_particle_force_on_fluid << std::endl;
  deallog << "Total particle force on the particles "
          << total_particle_force_on_particles << std::endl;

  if (output_vtu)
    {
      DataOut<3> data_out;
      data_out.add_data_vector(particle_projector.dof_handler,
                               particle_projector.void_fraction_solution,
                               "void_fraction");
      data_out.add_data_vector(
        particle_projector.fluid_force_on_particles_two_way_coupling
          .dof_handler,
        particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution,
        "force_pf");
      data_out.build_patches();

      data_out.write_vtu_with_pvtu_record(
        "./", vtu_label, number_quadrature_points, MPI_COMM_WORLD, 2, 8);
    }
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;
      deallog
        << "Particle-fluid Constant Force: fe_degree=1    number_quadrature_points=2"
        << std::endl;
      Functions::ConstantFunction<3> constant_func({0, 1., 0});
      test_void_fraction_qcm(1, 2, constant_func, false, "constant");

      deallog
        << "Particle-fluid Constant Force: fe_degree=1    number_quadrature_points=3"
        << std::endl;
      test_void_fraction_qcm(1, 3, constant_func, false, "constant");

      deallog
        << "Particle-fluid Linear Force: fe_degree=1    number_quadrature_points=2"
        << std::endl;
      LinearForce<3> linear_func;
      test_void_fraction_qcm(1, 2, linear_func, false, "linear");

      deallog
        << "Particle-fluid Linear Force: fe_degree=1    number_quadrature_points=3"
        << std::endl;
      test_void_fraction_qcm(1, 3, linear_func, false, "linear");
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
