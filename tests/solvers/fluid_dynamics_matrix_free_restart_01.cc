// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the matrix-free fluid dynamics solver restores
 * its solution vectors with up to date ghost values when it restarts, so that
 * the dynamic flow control of the first restarted time step calculates the
 * same beta force as an uninterrupted run. It follows the same stages as
 * cfd_dem_matrix_free_restart_01:
 *   1. carries out three time steps,
 *   2. writes a checkpoint,
 *   3. carries out three more time steps in the same solver, which gives the
 *      reference state of a continuous run,
 *   4. builds a second solver on the coarse mesh, whose fields are all zero,
 *   5. restarts that second solver from the checkpoint through
 *      set_initial_condition(), which is the production path,
 *   6. carries out the same three time steps as in stage 3.
 *
 * A signature of the state is printed after each of these stages. The
 * signature printed after stage 5 must match the one printed at stage 1, and
 * the one printed after stage 6 must match the one printed at stage 3. The
 * average velocity and the beta force calculated by the flow control are also
 * printed at every time step, and those of the restarted time steps must match
 * those of stage 3. Stage 4 is what gives the test its value: without a solver
 * whose state is empty to begin with, the values printed after the restart
 * would also be obtained by a restart that restores nothing at all.
 *
 * The signature is made of the continuous L2 norms of the present and previous
 * solutions, and of the average velocity at the inlet of the flow control.
 * These quantities are independent of the parallel partitioning, so the serial
 * and the mpirun=2 outputs of this test hold the same values. They are
 * evaluated over the locally owned cells and therefore read the ghost values of
 * the vectors, as the flow control does at the beginning of every time step. A
 * restart that leaves the ghost values of the restored vectors stale is
 * detected under mpirun=2: in Debug, deal.II stops on the assertion that
 * forbids reading a ghost entry of a vector that has not imported its ghost
 * values; in Release, the ghost entries read as zero, which lowers the norms
 * and the average velocity printed after the restart and changes the beta
 * force of the first restarted time step. With two processes, the mesh is
 * split at y = 0.5, a cut that crosses the inlet x = 0, so the inlet holds
 * degrees of freedom owned by the other process.
 *
 * The fluid solution is imposed analytically rather than solved, so that the
 * state of the simulation only depends on the time step at which it is
 * evaluated. This is what allows the restarted run to be compared to the
 * continuous one. The analytical solution is linear in space and the finite
 * element is of degree two, so the solution belongs to the finite element
 * space and the printed quantities do not depend on the mesh nor on the
 * partitioning. Its magnitude grows in time, so the present solution and the
 * previous solutions all differ from one another.
 *
 * The target average velocity of the flow control is larger than the average
 * velocity at the inlet at every time step of the test. The main controller is
 * therefore active from the third time step on and never clips the beta force
 * to zero, so that an error on the average velocity always changes the beta
 * force.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/parameters.h>

#include <solvers/fluid_dynamics_matrix_free.h>
#include <solvers/postprocessing_cfd.h>
#include <solvers/simulation_parameters.h>

// Tests
#include <../tests/tests.h>

// Standard library
#include <string>

/// Number of time steps carried out before the checkpoint and after it.
constexpr unsigned int number_of_steps = 3;

/// Boundary id of the inlet of the flow control: the face x = 0 of the
/// colorized hyper cube.
constexpr types::boundary_id inlet_boundary_id = 0;

/// Target average velocity of the flow control. It is larger than the average
/// velocity of the analytical solution at the inlet, 2 (1 + t / 2), at every
/// time step of the test.
constexpr double target_average_velocity = 4.;

/**
 * @brief Analytical solution imposed on the fluid. Its spatial shape is fixed
 * and linear and it is scaled by a time-dependent factor, as in
 * average_velocities_05. Its average velocity at the inlet x = 0 of the unit
 * square is twice the scaling factor.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class ScaledSolution : public Function<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] scale Factor that multiplies the spatial shape.
   */
  ScaledSolution(const double scale)
    : Function<dim>(dim + 1)
    , scale(scale)
  {}

  /**
   * @brief Evaluate the velocity and the pressure at a point.
   *
   * @param[in] p Point at which the solution is evaluated.
   * @param[out] values Velocity components followed by the pressure.
   */
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;

private:
  /// Factor that multiplies the spatial shape.
  const double scale;
};

template <int dim>
void
ScaledSolution<dim>::vector_value(const Point<dim> &p,
                                  Vector<double>   &values) const
{
  const double x = p[0];
  const double y = p[1];

  values(0) = scale * (1.0 + x + 2.0 * y);
  values(1) = scale * (2.0 - 3.0 * x + y);
  if constexpr (dim == 3)
    values(2) = scale * (0.25 + p[2] - x);

  values(dim) = scale * (0.5 + x - y);
}


/**
 * @brief Matrix-free fluid dynamics solver driven step by step by the test.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class FluidDynamicsMatrixFreeRestart : public FluidDynamicsMatrixFree<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] nsparam Parameters of the simulation.
   */
  FluidDynamicsMatrixFreeRestart(SimulationParameters<dim> &nsparam)
    : FluidDynamicsMatrixFree<dim>(nsparam)
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
  using VectorType = LinearAlgebra::distributed::Vector<double>;

  /**
   * @brief Set up the solver on the triangulation that has already been
   * created by the caller.
   */
  void
  setup();

  /**
   * @brief Fill the present solution with the analytical solution evaluated
   * at a given time.
   *
   * @param[in] time Time at which the analytical solution is evaluated.
   */
  void
  set_solution(const double time);

  /**
   * @brief Calculate the continuous L2 norm of a finite element field. It
   * reads the ghost values of the field and does not depend on the mesh nor on
   * the parallel partitioning.
   *
   * @param[in] field Finite element field of the fluid dynamics.
   *
   * @return Continuous L2 norm of the field.
   */
  double
  field_l2_norm(const VectorType &field);

  /**
   * @brief Calculate the average velocity of the present solution at the
   * inlet with the function used by the dynamic flow control.
   *
   * @return Average velocity at the inlet.
   */
  double
  inlet_average_velocity();

  /**
   * @brief Print the signature of the state of the simulation. The ghost
   * values of the vectors are deliberately not refreshed here: a vector whose
   * ghost values are left stale by the restart must be detected rather than
   * repaired.
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

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::setup()
{
  this->setup_dofs();

  // The distributed vectors of the matrix-free solver refuse to have a ghost
  // entry read while they are not in their ghosted state. The ones that have
  // just been initialized are all zero and have never been ghosted, so their
  // ghost values are imported once here to make the state of the empty solver
  // printable.
  this->present_solution->update_ghost_values();
  for (auto &previous_solution : *this->previous_solutions)
    previous_solution.update_ghost_values();
}

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::set_solution(const double time)
{
  const ScaledSolution<dim> solution(1. + 0.5 * time);

  // VectorTools::interpolate writes the degrees of freedom of every locally
  // owned cell, including those owned by another process, so it needs a vector
  // with ghost entries. As set_nodal_values() does, the interpolation is
  // carried out in the newton update, which is built from the partitioner of
  // the matrix-free operator, and the result is copied into the present
  // solution, whose ghost values are then refreshed as they are after a solve.
  VectorTools::interpolate(*this->get_mapping(),
                           *this->dof_handler,
                           solution,
                           this->newton_update);

  *this->present_solution = this->newton_update;
  this->present_solution->update_ghost_values();
}

template <int dim>
double
FluidDynamicsMatrixFreeRestart<dim>::field_l2_norm(const VectorType &field)
{
  Vector<double> cellwise_norm(this->triangulation->n_active_cells());

  VectorTools::integrate_difference(*this->get_mapping(),
                                    *this->dof_handler,
                                    field,
                                    Functions::ZeroFunction<dim>(dim + 1),
                                    cellwise_norm,
                                    *this->cell_quadrature,
                                    VectorTools::L2_norm);

  return VectorTools::compute_global_error(*this->triangulation,
                                           cellwise_norm,
                                           VectorTools::L2_norm);
}

template <int dim>
double
FluidDynamicsMatrixFreeRestart<dim>::inlet_average_velocity()
{
  return calculate_average_velocity(*this->dof_handler,
                                    *this->present_solution,
                                    inlet_boundary_id,
                                    *this->face_quadrature,
                                    *this->get_mapping());
}

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::print_state(const std::string &label)
{
  deallog << label << std::endl;
  deallog << "  Time                   : "
          << this->simulation_control->get_current_time() << std::endl;
  deallog << "  Iteration              : "
          << this->simulation_control->get_iteration_number() << std::endl;
  deallog << "  Present solution       : "
          << field_l2_norm(*this->present_solution) << std::endl;
  for (unsigned int i = 0; i < this->previous_solutions->size(); ++i)
    deallog << "  Previous solution " << i << "    : "
            << field_l2_norm((*this->previous_solutions)[i]) << std::endl;
  deallog << "  Inlet average velocity : " << inlet_average_velocity()
          << std::endl;
}

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::advance(
  const unsigned int number_of_time_steps)
{
  const unsigned int flow_direction =
    this->simulation_parameters.flow_control.flow_direction;

  for (unsigned int step = 0; step < number_of_time_steps; ++step)
    {
      this->simulation_control->integrate();

      // This is the production function, called at the beginning of the time
      // step as in FluidDynamicsMatrixFree::solve(). It calculates the beta
      // force from the average velocity of the present solution, which is
      // still the solution of the previous time step, the restored one on the
      // first time step after a restart.
      this->dynamic_flow_control();

      deallog << "Time step "
              << this->simulation_control->get_iteration_number()
              << " : average velocity = " << inlet_average_velocity()
              << ", beta = " << this->flow_control.get_beta()[flow_direction]
              << std::endl;

      set_solution(this->simulation_control->get_current_time());
      this->finish_time_step();
    }
}

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::run_and_checkpoint()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1, true);
  this->triangulation->refine_global(3);

  setup();

  // Initial condition, set as the nodal initial condition of
  // FluidDynamicsMatrixFree::set_initial_condition_fd() is: the solution is
  // imposed and finish_time_step() copies it into the previous solutions.
  set_solution(this->simulation_control->get_current_time());
  this->finish_time_step();

  advance(number_of_steps);
  print_state("State before the checkpoint");

  this->write_checkpoint();

  // Keep going in the same solver. This is the reference that the restarted
  // run must reproduce.
  advance(number_of_steps);
  print_state("State at the end of the continuous run");
}

template <int dim>
void
FluidDynamicsMatrixFreeRestart<dim>::restart_and_check()
{
  // Only the coarse mesh is created here. The refinement of the triangulation
  // comes from the checkpoint, which is also what read_mesh_and_manifolds does
  // on a restart since it skips the initial refinement in that case.
  GridGenerator::hyper_cube(*this->triangulation, 0, 1, true);

  setup();
  print_state("State before the restart");

  // Restart the simulation through the production path.
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type, true);
  print_state("State after the restart");

  advance(number_of_steps);
  print_state("State at the end of the restarted run");
}

/**
 * @brief Set the parameters of the solvers of the test.
 *
 * @tparam dim Number of spatial dimensions.
 *
 * @param[out] parameters Parameters of the simulation.
 * @param[in] restart Whether the simulation is restarted from a checkpoint.
 */
template <int dim>
void
setup_parameters(SimulationParameters<dim> &parameters, const bool restart)
{
  // Every boundary id of the triangulation needs a boundary condition. The
  // parsed default boundary conditions are of type none on the boundary ids 0
  // to 2 dim - 1 of the colorized hyper cube. The boundary id 0 is the inlet of
  // the flow control.
  ParameterHandler              prm;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 2 * dim;
  size_of_subsections.manifolds           = 0;

  parameters.declare(prm, size_of_subsections);
  parameters.parse(prm);

  parameters.simulation_control.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf2;
  parameters.simulation_control.dt = 0.1;
  parameters.simulation_control.end_control =
    Parameters::SimulationControl::EndControl::time;
  parameters.simulation_control.time_end = 1.0;

  // A degree two element is used so that the imposed solution, which is
  // linear in space, belongs to the finite element space. The matrix-free
  // solver requires an equal order interpolation.
  parameters.fem_parameters.velocity_degree = 2;
  parameters.fem_parameters.pressure_degree = 2;

  // The checkpointing of finish_time_step() is disabled: the test calls
  // write_checkpoint() itself at the time step of its choosing, which avoids
  // overwriting the checkpoint with the time steps that follow it.
  parameters.restart_parameters.checkpoint = false;
  parameters.restart_parameters.restart    = restart;
  parameters.restart_parameters.frequency  = 1;
  parameters.restart_parameters.filename =
    "fluid_dynamics_matrix_free_restart_01_restart";

  // A threshold of zero keeps every change of the beta force.
  parameters.flow_control.enable_flow_control = true;
  parameters.flow_control.boundary_flow_id    = inlet_boundary_id;
  parameters.flow_control.flow_direction      = 0;
  parameters.flow_control.average_velocity_0  = target_average_velocity;
  parameters.flow_control.beta_0              = 1.;
  parameters.flow_control.alpha               = 1.;
  parameters.flow_control.beta_threshold      = 0.;
  parameters.flow_control.verbosity           = Parameters::Verbosity::quiet;

  parameters.mesh_adaptation.type = Parameters::MeshAdaptation::Type::none;
  parameters.timer.type           = Parameters::Timer::Type::none;

  parameters.linear_solver.at(PhysicsID::fluid_dynamics).verbosity =
    Parameters::Verbosity::quiet;
  parameters.physics_solving_strategy.at(PhysicsID::fluid_dynamics).verbosity =
    Parameters::Verbosity::quiet;

  Parameters::PhysicalProperties physical_properties;
  physical_properties.fluids.resize(1);
  physical_properties.number_of_fluids = 1;
  physical_properties.fluids[0].rheological_model =
    Parameters::Material::RheologicalModel::newtonian;
  physical_properties.fluids[0].kinematic_viscosity = 1;
  physical_properties.fluids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.number_of_solids                = 0;
  physical_properties.number_of_material_interactions = 0;

  parameters.physical_properties_manager.initialize(physical_properties);
}

void
test()
{
  {
    SimulationParameters<2> parameters;
    setup_parameters(parameters, false);

    FluidDynamicsMatrixFreeRestart<2> continuous_run(parameters);
    continuous_run.run_and_checkpoint();
  }

  // Make sure that every checkpoint file has been written before any process
  // starts reading them back.
  MPI_Barrier(MPI_COMM_WORLD);

  {
    SimulationParameters<2> parameters;
    setup_parameters(parameters, true);

    FluidDynamicsMatrixFreeRestart<2> restarted_run(parameters);
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
