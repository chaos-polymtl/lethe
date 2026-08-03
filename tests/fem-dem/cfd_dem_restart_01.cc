// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the state of the matrix-based unresolved CFD-DEM
 * solver is restored by a restart. The three families of fields that the
 * checkpoint of this solver carries are monitored: the fluid solution and its
 * BDF history, the void fraction and its BDF history, and the particles.
 *
 * The test:
 *   1. carries out three time steps,
 *   2. writes a checkpoint,
 *   3. carries out three more time steps in the same solver, which gives the
 *      reference state of a continuous run,
 *   4. builds a second solver on the coarse mesh, which has no particles and
 *      whose fields are all zero,
 *   5. restarts that second solver from the checkpoint,
 *   6. carries out the same three time steps as in stage 3.
 *
 * A signature of the state is printed after each of these stages. The signature
 * printed after stage 5 must match the one printed at stage 1, and the one
 * printed after stage 6 must match the one printed at stage 3.
 *
 * Stage 4 is what gives the test its value: without a solver whose state is
 * empty to begin with, the values printed after the restart would also be
 * obtained by a restart that restores nothing at all.
 *
 * The restart is carried out with read_checkpoint(), which is the production
 * path, so the fields go through the same SolutionTransfer::deserialize calls,
 * the same triangulation load and the same ParticleHandler::deserialize as in a
 * simulation.
 *
 * The fluid solution is imposed analytically and the particles are moved with
 * their own velocity rather than by the DEM force loop, so that the state of
 * the simulation only depends on the time step at which it is evaluated. This
 * is what allows the restarted run to be compared to the continuous one. Only
 * the void fraction is calculated by the solver itself, from the position of
 * the particles.
 *
 * The particles all follow different trajectories and the void fraction is
 * calculated with the quadrature centered method, whose result varies
 * continuously with the position of the particles. The fields of the two time
 * histories therefore all differ from one another, which is what makes the
 * test able to tell a history that is restored correctly from one whose fields
 * are mixed up, dropped or duplicated.
 *
 * The signature is made of continuous L2 norms and of global particle
 * checksums, both of which are independent of the parallel partitioning, so the
 * serial and the mpirun=2 outputs of this test are identical. Since the L2
 * norms are assembled over the locally owned cells, they read the ghost values
 * of the fields, and a field whose ghost values are left stale by the restart
 * is detected under mpirun=2.
 *
 * The present solution and the first previous solution carry the same values at
 * the time of the checkpoint. This is not a weakness of the test but a property
 * of the time loop of every Lethe solver: the checkpoint is written by
 * finish_time_step(), right after the time vectors have been percolated. The
 * second previous solution differs from both, so the transfer of the time
 * history is covered.
 */

// Deal.II includes
#include <deal.II/base/parameter_handler.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/dem_properties.h>

#include <fem-dem/cfd_dem_coupling.h>

// Tests
#include <../tests/fem-dem/cfd_dem_test_utilities.h>

#include <../tests/tests.h>

/// Diameter of the particles of the test. It is large enough for the solid
/// fraction of the cells that hold a particle to be of the order of ten
/// percent, so that the void fraction field is far from uniform.
constexpr double particle_diameter = 0.15;

/// Density of the particles of the test. It is only used to establish their
/// mass, which is checkpointed but does not influence the void fraction.
constexpr double particle_density = 1000.;

/// Number of time steps carried out before the checkpoint and after it.
constexpr unsigned int number_of_steps = 3;


template <int dim, typename PropertiesIndex>
class CFDDEMRestart : public CFDDEMSolver<dim, PropertiesIndex>
{
public:
  CFDDEMRestart(CFDDEMSimulationParameters<dim> &nsparam)
    : CFDDEMSolver<dim, PropertiesIndex>(nsparam)
  {}

  /**
   * @brief Carry out the time steps of the continuous run and write the
   * checkpoint that the restarted run reads.
   */
  void
  run_and_checkpoint();

  /**
   * @brief Restart from that checkpoint and carry out the same time steps as
   * the end of the continuous run.
   */
  void
  restart_and_check();

private:
  /**
   * @brief Set up the solver on the triangulation that has already been
   * created by the caller. This is the sequence of calls that
   * CFDDEMSolver::solve() carries out between the reading of the mesh and the
   * setting of the initial condition.
   */
  void
  setup();

  /**
   * @brief Fill the fluid solution with the analytical solution evaluated at a
   * given time.
   *
   * @param[in] time Time at which the analytical solution is evaluated.
   */
  void
  set_fluid_solution(const double time);

  /**
   * @brief Print the signature of the state of the simulation.
   *
   * @param[in] label Description of the point of the simulation at which the
   * signature is printed.
   */
  void
  print_state(const std::string &label);

  /**
   * @brief Carry out a given number of time steps.
   *
   * @param[in] number_of_time_steps Number of time steps to carry out.
   */
  void
  advance(const unsigned int number_of_time_steps);
};

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::setup()
{
  this->dem_setup_parameters();
  this->setup_dofs();
  this->vertices_cell_mapping();
}

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::set_fluid_solution(const double time)
{
  const ScaledFluidSolution<dim> solution(1. + 0.5 * time);

  // VectorTools::interpolate writes the locally owned entries, so the
  // interpolation is carried out in a vector without ghost entries and the
  // result is copied into the vectors of the solver, which carry the ghost
  // values expected by the postprocessing.
  GlobalVectorType interpolated(this->locally_owned_dofs,
                                this->mpi_communicator);
  VectorTools::interpolate(*this->get_mapping(),
                           *this->dof_handler,
                           solution,
                           interpolated);

  this->local_evaluation_point = interpolated;
  *this->present_solution      = interpolated;
}

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::print_state(const std::string &label)
{
  deallog << label << std::endl;
  deallog << "  Time                      : "
          << this->simulation_control->get_current_time() << std::endl;
  deallog << "  Iteration                 : "
          << this->simulation_control->get_iteration_number() << std::endl;

  deallog << "  Fluid solution            : "
          << field_l2_norm(*this->get_mapping(),
                           *this->dof_handler,
                           *this->cell_quadrature,
                           *this->present_solution,
                           dim + 1)
          << std::endl;
  for (unsigned int i = 0; i < this->previous_solutions->size(); ++i)
    deallog << "  Previous fluid solution " << i << " : "
            << field_l2_norm(*this->get_mapping(),
                             *this->dof_handler,
                             *this->cell_quadrature,
                             (*this->previous_solutions)[i],
                             dim + 1)
            << std::endl;

  deallog
    << "  Void fraction             : "
    << field_l2_norm(*this->particle_projector.mapping,
                     this->particle_projector.dof_handler,
                     *this->particle_projector.quadrature,
                     this->particle_projector.void_fraction_locally_relevant,
                     1)
    << std::endl;
  for (unsigned int i = 0;
       i < this->particle_projector.previous_void_fraction.size();
       ++i)
    deallog << "  Previous void fraction " << i << "  : "
            << field_l2_norm(*this->particle_projector.mapping,
                             this->particle_projector.dof_handler,
                             *this->particle_projector.quadrature,
                             this->particle_projector.previous_void_fraction[i],
                             1)
            << std::endl;

  print_particle_state<dim, PropertiesIndex>(this->particle_handler,
                                             this->mpi_communicator);
}

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::advance(
  const unsigned int number_of_time_steps)
{
  for (unsigned int step = 0; step < number_of_time_steps; ++step)
    {
      this->simulation_control->integrate();

      const double time = this->simulation_control->get_current_time();

      set_fluid_solution(time);

      move_test_particles<dim, PropertiesIndex>(
        this->particle_handler, this->simulation_control->get_time_step());
      this->sort_particles_into_subdomains_and_cells();

      // calculate_void_fraction and finish_time_step_fd are the production
      // functions, which guarantees that the void fraction and the time
      // histories are established exactly as they are during a simulation.
      this->calculate_void_fraction(time);
      this->finish_time_step_fd();
    }
}

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::run_and_checkpoint()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);
  this->triangulation->refine_global(2);

  setup();

  // The particles are generated inside a box that is strictly contained in the
  // domain so that they remain in it for the whole duration of the test.
  Point<dim> bottom_left;
  Point<dim> top_right;
  for (int d = 0; d < dim; ++d)
    {
      bottom_left[d] = 0.15;
      top_right[d]   = 0.85;
    }

  insert_test_particles<dim, PropertiesIndex>(
    dynamic_cast<parallel::distributed::Triangulation<dim> &>(
      *this->triangulation),
    this->particle_handler,
    bottom_left,
    top_right,
    1,
    particle_diameter,
    particle_density);

  this->sort_particles_into_subdomains_and_cells();
  this->particle_projector.initialize_void_fraction(
    this->simulation_control->get_current_time());

  advance(number_of_steps);
  print_state("State before the checkpoint");

  this->write_checkpoint();

  // Keep going in the same solver. This is the reference that the restarted
  // run must reproduce.
  advance(number_of_steps);
  print_state("State at the end of the continuous run");
}

template <int dim, typename PropertiesIndex>
void
CFDDEMRestart<dim, PropertiesIndex>::restart_and_check()
{
  // Only the coarse mesh is created here. The refinement of the triangulation
  // comes from the checkpoint, which is also what read_mesh_and_manifolds does
  // on a restart since it skips the initial refinement in that case.
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);

  setup();
  print_state("State before the restart");

  this->read_checkpoint();

  // The vertex to cell map and the cells of the particles are rebuilt on the
  // mesh that has just been read, as initialize_dem_parameters() does after the
  // restart of a simulation.
  this->vertices_cell_mapping();
  this->sort_particles_into_subdomains_and_cells();

  print_state("State after the restart");

  advance(number_of_steps);
  print_state("State at the end of the restarted run");
}

void
test()
{
  using PropertiesIndex = DEM::CFDDEMProperties::PropertiesIndex;

  const std::string restart_filename = "cfd_dem_restart_01_restart";

  {
    CFDDEMSimulationParameters<3> parameters;
    setup_cfd_dem_test_parameters<3>(parameters,
                                     particle_diameter,
                                     restart_filename,
                                     false);

    CFDDEMRestart<3, PropertiesIndex> continuous_run(parameters);
    continuous_run.run_and_checkpoint();
  }

  // Make sure that every checkpoint file has been written before any process
  // starts reading them back.
  MPI_Barrier(MPI_COMM_WORLD);

  {
    CFDDEMSimulationParameters<3> parameters;
    setup_cfd_dem_test_parameters<3>(parameters,
                                     particle_diameter,
                                     restart_filename,
                                     true);

    CFDDEMRestart<3, PropertiesIndex> restarted_run(parameters);
    restarted_run.restart_and_check();
  }
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
