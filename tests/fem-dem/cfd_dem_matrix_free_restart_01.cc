// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the state of the matrix-free unresolved CFD-DEM
 * solver is restored by a restart. It is the matrix-free counterpart of
 * cfd_dem_restart_01 and follows exactly the same stages:
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
 * printed after stage 6 must match the one printed at stage 3. Stage 4 is what
 * gives the test its value: without a solver whose state is empty to begin
 * with, the values printed after the restart would also be obtained by a
 * restart that restores nothing at all.
 *
 * The void fraction of this test is monitored through the deal.II distributed
 * vectors of the projector, which are the ones the matrix-free operator
 * consumes, rather than through the Trilinos vectors monitored by
 * cfd_dem_restart_01. These are the vectors the checkpoint of this solver
 * carries and, contrary to the Trilinos ones, they are restored in the same way
 * whether or not the library is built with LETHE_USE_LDV.
 *
 * The time derivative of the void fraction is evaluated at every time step, as
 * the production time loop does, so that the matrix-free consumer of the void
 * fraction history is exercised. It is deliberately left out of the printed
 * signature: it is not checkpointed but reevaluated from the void fraction and
 * its history at the beginning of every step, so it is legitimately zero right
 * after a restart and comparing it there would be comparing a quantity that the
 * checkpoint never carried.
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

#include <fem-dem/cfd_dem_coupling_matrix_free.h>

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
class CFDDEMMatrixFreeRestart : public CFDDEMMatrixFree<dim, PropertiesIndex>
{
public:
  CFDDEMMatrixFreeRestart(CFDDEMSimulationParameters<dim> &nsparam)
    : CFDDEMMatrixFree<dim, PropertiesIndex>(nsparam)
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
   * CFDDEMMatrixFree::solve() carries out between the reading of the mesh and
   * the setting of the initial condition. Contrary to the matrix-based solver,
   * the matrix-free solver has no vertex to cell map of its own: the one used
   * by the quadrature centered method belongs to the projector and is built by
   * its setup_dofs().
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
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::setup()
{
  this->dem_setup_parameters();
  this->setup_dofs();

  // The distributed vectors of the matrix-free solver refuse to have a ghost
  // entry read while they are not in their ghosted state. The ones that have
  // just been initialized are all zero and have never been ghosted, so their
  // ghost values are imported once here to make the state of the empty solver
  // printable. print_state() deliberately does not do this itself: a field
  // whose ghost values are left stale by the restart must be detected rather
  // than repaired.
  this->present_solution->update_ghost_values();
  for (auto &previous_solution : *this->previous_solutions)
    previous_solution.update_ghost_values();

  this->particle_projector.void_fraction_solution.update_ghost_values();
  for (auto &previous_void_fraction :
       this->particle_projector.void_fraction_previous_solution)
    previous_void_fraction.update_ghost_values();
}

template <int dim, typename PropertiesIndex>
void
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::set_fluid_solution(
  const double time)
{
  const ScaledFluidSolution<dim> solution(1. + 0.5 * time);

  // The vectors of the matrix-free solver are built from the partitioner of
  // its operator and therefore carry ghost entries. VectorTools::interpolate
  // only writes the locally owned entries, so the interpolation is carried out
  // in the newton update, as set_nodal_values() does, and the ghost values of
  // the present solution are refreshed afterwards.
  VectorTools::interpolate(*this->get_mapping(),
                           *this->dof_handler,
                           solution,
                           this->newton_update);

  this->local_evaluation_point = this->newton_update;
  *this->present_solution      = this->newton_update;
  this->present_solution->update_ghost_values();
}

template <int dim, typename PropertiesIndex>
void
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::print_state(
  const std::string &label)
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

  deallog << "  Void fraction             : "
          << field_l2_norm(*this->particle_projector.mapping,
                           this->particle_projector.dof_handler,
                           *this->particle_projector.quadrature,
                           this->particle_projector.void_fraction_solution,
                           1)
          << std::endl;
  for (unsigned int i = 0;
       i < this->particle_projector.void_fraction_previous_solution.size();
       ++i)
    deallog << "  Previous void fraction " << i << "  : "
            << field_l2_norm(
                 *this->particle_projector.mapping,
                 this->particle_projector.dof_handler,
                 *this->particle_projector.quadrature,
                 this->particle_projector.void_fraction_previous_solution[i],
                 1)
            << std::endl;

  print_particle_state<dim, PropertiesIndex>(this->particle_handler,
                                             this->mpi_communicator);
}

template <int dim, typename PropertiesIndex>
void
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::advance(
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

      // These are the production functions, which guarantees that the void
      // fraction, its time derivative and the time histories are established
      // exactly as they are during a simulation. Contrary to the matrix-based
      // solver, the void fraction of the matrix-free solver is calculated by
      // the projector itself rather than through a function of the solver.
      this->particle_projector.calculate_void_fraction(time);
      this->evaluate_time_derivative_void_fraction();
      this->finish_time_step_fd();
    }
}

template <int dim, typename PropertiesIndex>
void
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::run_and_checkpoint()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);
  this->triangulation->refine_global(2);

  setup();

  // The particles are generated inside a box that is strictly contained in the
  // domain so that they remain in it for the whole duration of the test.
  Point<dim> bottom_left;
  Point<dim> top_right;
  for (unsigned int d = 0; d < dim; ++d)
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
CFDDEMMatrixFreeRestart<dim, PropertiesIndex>::restart_and_check()
{
  // Only the coarse mesh is created here. The refinement of the triangulation
  // comes from the checkpoint, which is also what read_mesh_and_manifolds does
  // on a restart since it skips the initial refinement in that case.
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);

  setup();
  print_state("State before the restart");

  this->read_checkpoint();

  // The cells of the particles are rebuilt on the mesh that has just been
  // read, as initialize_dem_parameters() does after the restart of a
  // simulation.
  this->sort_particles_into_subdomains_and_cells();

  print_state("State after the restart");

  advance(number_of_steps);
  print_state("State at the end of the restarted run");
}

void
test()
{
  using PropertiesIndex = DEM::CFDDEMProperties::PropertiesIndex;

  const std::string restart_filename = "cfd_dem_matrix_free_restart_01_restart";

  {
    CFDDEMSimulationParameters<3> parameters;
    setup_cfd_dem_test_parameters<3>(parameters,
                                     particle_diameter,
                                     restart_filename,
                                     false);

    CFDDEMMatrixFreeRestart<3, PropertiesIndex> continuous_run(parameters);
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

    CFDDEMMatrixFreeRestart<3, PropertiesIndex> restarted_run(parameters);
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
