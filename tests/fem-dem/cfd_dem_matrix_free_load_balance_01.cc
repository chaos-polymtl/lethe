// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the state of the matrix-free unresolved CFD-DEM
 * solver survives a load balance. It is the matrix-free counterpart of
 * cfd_dem_load_balance_01 and follows exactly the same stages:
 *   1. carries out three time steps,
 *   2. repartitions the triangulation with load_balance(),
 *   3. carries out three more time steps.
 *
 * A signature of the state is printed after each of these stages. The
 * signatures printed before and after the load balance must be identical: the
 * repartitioning moves cells and particles between the processes but changes
 * neither the mesh nor the finite element space, and the imposed fields belong
 * to that space.
 *
 * The load balance is carried out with load_balance(), which is the production
 * path, so the fields go through the same
 * SolutionTransfer::prepare_for_coarsening_and_refinement and
 * ParticleHandler::prepare_for_coarsening_and_refinement calls as in a
 * simulation. It is triggered at every call by the load balancing parameters
 * set in cfd_dem_test_utilities.h, so the test does not depend on the iteration
 * at which it asks for it.
 *
 * The void fraction of this test is monitored through the deal.II distributed
 * vectors of the projector, which are the ones the matrix-free operator
 * consumes, rather than through the Trilinos vectors monitored by
 * cfd_dem_load_balance_01. These are the vectors the load balance of this
 * solver transfers and, contrary to the Trilinos ones, they are restored in the
 * same way whether or not the library is built with LETHE_USE_LDV.
 *
 * The time derivative of the void fraction is evaluated at every time step, as
 * the production time loop does, so that the matrix-free consumer of the void
 * fraction history is exercised. It is deliberately left out of the printed
 * signature: it is not transferred but reevaluated from the void fraction and
 * its history at the beginning of every step, so it is legitimately left
 * untouched by the repartitioning.
 *
 * The fluid solution is imposed analytically and the particles are moved with
 * their own velocity rather than by the DEM force loop, so that the state of
 * the simulation only depends on the time step at which it is evaluated. Only
 * the void fraction is calculated by the solver itself, from the position of
 * the particles.
 *
 * The particles all follow different trajectories and the void fraction is
 * calculated with the quadrature centered method, whose result varies
 * continuously with the position of the particles. The fields of the time
 * history therefore all differ from one another, which is what makes the test
 * able to tell a history that is transferred correctly from one whose fields
 * are mixed up, dropped or duplicated.
 *
 * The signature is made of continuous L2 norms and of global particle
 * checksums, both of which are independent of the parallel partitioning. This
 * is what makes them comparable across a repartitioning, and it also makes the
 * serial and the mpirun=2 outputs of this test identical. In serial the
 * repartitioning has nothing to move, but the packing and the interpolation of
 * every field are still exercised.
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
/// mass, which is transferred but does not influence the void fraction.
constexpr double particle_density = 1000.;

/// Number of time steps carried out before the load balance and after it.
constexpr unsigned int number_of_steps = 3;


template <int dim, typename PropertiesIndex>
class CFDDEMMatrixFreeLoadBalance
  : public CFDDEMMatrixFree<dim, PropertiesIndex>
{
public:
  CFDDEMMatrixFreeLoadBalance(CFDDEMSimulationParameters<dim> &nsparam)
    : CFDDEMMatrixFree<dim, PropertiesIndex>(nsparam)
  {}

  /**
   * @brief Carry out the time steps of the simulation and repartition the
   * triangulation in the middle of them.
   */
  void
  run();

private:
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
CFDDEMMatrixFreeLoadBalance<dim, PropertiesIndex>::set_fluid_solution(
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
CFDDEMMatrixFreeLoadBalance<dim, PropertiesIndex>::print_state(
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
CFDDEMMatrixFreeLoadBalance<dim, PropertiesIndex>::advance(
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
CFDDEMMatrixFreeLoadBalance<dim, PropertiesIndex>::run()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);
  this->triangulation->refine_global(2);

  // Contrary to the matrix-based solver, the matrix-free solver has no vertex
  // to cell map of its own: the one used by the quadrature centered method
  // belongs to the projector and is built by its setup_dofs().
  this->dem_setup_parameters();
  this->setup_dofs();

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

  // initialize_dem_parameters sorts the particles into their cells and builds
  // the cell neighbors and the boundary cells that load_balance() updates.
  this->initialize_dem_parameters();
  this->particle_projector.initialize_void_fraction(
    this->simulation_control->get_current_time());

  advance(number_of_steps);
  print_state("State before the load balance");

  this->load_balance();
  print_state("State after the load balance");

  // Keep going on the new partition. If the transfer had lost a field, the
  // simulation would carry on from a state that differs from the one it had
  // before the repartitioning.
  advance(number_of_steps);
  print_state("State at the end of the simulation");
}

void
test()
{
  using PropertiesIndex = DEM::CFDDEMProperties::PropertiesIndex;

  const std::string restart_filename =
    "cfd_dem_matrix_free_load_balance_01_restart";

  CFDDEMSimulationParameters<3> parameters;
  setup_cfd_dem_test_parameters<3>(parameters,
                                   particle_diameter,
                                   restart_filename,
                                   false);

  CFDDEMMatrixFreeLoadBalance<3, PropertiesIndex> problem(parameters);
  problem.run();
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
