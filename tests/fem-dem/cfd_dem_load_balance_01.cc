// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the state of the matrix-based unresolved CFD-DEM
 * solver survives a load balance. The same three families of fields as in
 * cfd_dem_restart_01 are monitored: the fluid solution and its BDF history, the
 * void fraction and its BDF history, and the particles.
 *
 * The test:
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
 * The fluid solution is imposed analytically and the particles are moved with
 * their own velocity rather than by the DEM force loop, so that the state of
 * the simulation only depends on the time step at which it is evaluated. Only
 * the void fraction is calculated by the solver itself, from the position of
 * the particles.
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

#include <fem-dem/cfd_dem_coupling.h>

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
class CFDDEMLoadBalance : public CFDDEMSolver<dim, PropertiesIndex>
{
public:
  CFDDEMLoadBalance(CFDDEMSimulationParameters<dim> &nsparam)
    : CFDDEMSolver<dim, PropertiesIndex>(nsparam)
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
CFDDEMLoadBalance<dim, PropertiesIndex>::set_fluid_solution(const double time)
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
CFDDEMLoadBalance<dim, PropertiesIndex>::print_state(const std::string &label)
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
CFDDEMLoadBalance<dim, PropertiesIndex>::advance(
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
CFDDEMLoadBalance<dim, PropertiesIndex>::run()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);
  this->triangulation->refine_global(2);

  this->dem_setup_parameters();
  this->setup_dofs();

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

  // initialize_dem_parameters sorts the particles into their cells and builds
  // the cell neighbors and the boundary cells that load_balance() updates.
  this->initialize_dem_parameters();
  this->vertices_cell_mapping();
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

  CFDDEMSimulationParameters<3> parameters;
  setup_cfd_dem_test_parameters<3>(parameters,
                                   particle_diameter,
                                   "cfd_dem_load_balance_01_restart",
                                   false);

  CFDDEMLoadBalance<3, PropertiesIndex> problem(parameters);
  problem.run();
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;

      // The L2 projection of the void fraction is solved with a SolverControl
      // that logs its history, and the residuals of its iterations depend on
      // the partitioning. Only the messages of the test itself, which carry the
      // "DEAL" and the MPI rank prefixes, are kept in the output file.
      deallog.depth_file(2);

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
