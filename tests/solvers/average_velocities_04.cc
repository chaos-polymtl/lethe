// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests that the time-averaged velocity and the Reynolds
 * stresses survive a mesh adaptation with deal.II vectors.
 *
 * The test derives from the matrix-free solver and:
 *   1. accumulates the averaged fields over five time steps,
 *   2. uniformly refines the mesh,
 *   3. accumulates the averaged fields over five more time steps.
 *
 * The continuous L2 norm of the three averaged fields is printed after each of
 * these stages. The norms printed before and after the refinement must be
 * identical.
 *
 * The refinement is carried out with refine_mesh_uniform(), which is the
 * the dynamic mesh adaptation of the solvers goes through, so
 * the averaged fields are transferred by the same calls to
 * AverageVelocities::prepare_for_mesh_adaptation and
 * AverageVelocities::post_mesh_adaptation as in a simulation.
 *
 * The monitored quantity is the continuous L2 norm rather than the l2_norm() of
 * the underlying vector because the former is independent of the mesh and of
 * the parallel partitioning. A uniform refinement embeds the coarse finite
 * element space into the refined one, so the transferred fields are identical
 * as functions and their L2 norms are unchanged. A failure of the transfer,
 * such as the reinitialization of the accumulators by setup_dofs() without
 * transferring them, makes those norms drop to zero.
 *
 * The analytical solution imposed at every time step is linear in space and the
 * finite element is of second order. The average velocity is then linear and
 * the Reynolds stresses, which are products of the velocity fluctuations, are
 * quadratic, so all the averaged fields belong to the finite element space and
 * their interpolants do not depend on the mesh. The norms printed at the end of
 * this test are therefore expected to be identical to the ones printed at the
 * end of average_velocities_05, which carries out the same time steps but
 * interrupts them with a restart instead of a mesh refinement.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/parameters.h>

#include <solvers/fluid_dynamics_matrix_free.h>
#include <solvers/simulation_parameters.h>

// Tests
#include <../tests/tests.h>

/**
 * @brief Analytical solution used to fill the solution vectors. The spatial
 * shape is fixed and is scaled by a time-dependent factor, so that the time
 * average differs from the instantaneous field and the Reynolds stresses do
 * not vanish.
 *
 * The spatial shape is linear so that, together with the second order finite
 * element used by the test, the averaged fields belong to the finite element
 * space. See the description at the top of the file.
 */
template <int dim>
class ScaledSolution : public Function<dim>
{
public:
  ScaledSolution(const double scale)
    : Function<dim>(dim + 1)
    , scale(scale)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;

private:
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


template <int dim>
class AverageVelocitiesMeshAdaptation : public FluidDynamicsMatrixFree<dim>
{
public:
  AverageVelocitiesMeshAdaptation(SimulationParameters<dim> nsparam)
    : FluidDynamicsMatrixFree<dim>(nsparam)
  {}

  void
  run();

private:
  using VectorType = LinearAlgebra::distributed::Vector<double>;

  /**
   * @brief Fill the solution vectors used by the averaging with the analytical
   * solution evaluated at a given time.
   *
   * @param[in] time Time at which the analytical solution is evaluated.
   */
  void
  set_solution(const double time);

  /**
   * @brief Calculate the continuous L2 norm of a finite element field. Contrary
   * to the l2_norm() of the underlying vector, this quantity does not depend on
   * the mesh nor on the parallel partitioning, which makes it comparable before
   * and after a mesh adaptation and between runs using a different number of
   * processes.
   *
   * @param[in] field Finite element field with up to date ghost values.
   *
   * @return Continuous L2 norm of the field.
   */
  double
  field_l2_norm(const VectorType &field);

  /**
   * @brief Print the L2 norm of the three averaged fields.
   *
   * @param[in] label Description of the point of the simulation at which the
   * norms are printed.
   */
  void
  print_averaged_fields(const std::string &label);

  /**
   * @brief Accumulate the averaged fields over a given number of time steps.
   *
   * @param[in] number_of_steps Number of time steps to carry out.
   */
  void
  advance(const unsigned int number_of_steps);
};

template <int dim>
void
AverageVelocitiesMeshAdaptation<dim>::set_solution(const double time)
{
  const ScaledSolution<dim> solution(1. + 0.5 * time);

  // VectorTools::interpolate writes the locally owned entries, so the
  // interpolation is carried out in a vector without ghost entries and the
  // result is copied into the vectors of the solver, which carry the ghost
  // values expected by the postprocessing.
  VectorType interpolated(this->locally_owned_dofs, this->mpi_communicator);
  VectorTools::interpolate(*this->get_mapping(),
                           *this->dof_handler,
                           solution,
                           interpolated);

  this->local_evaluation_point = interpolated;
  *this->present_solution      = interpolated;
  this->present_solution->update_ghost_values();
}

template <int dim>
double
AverageVelocitiesMeshAdaptation<dim>::field_l2_norm(const VectorType &field)
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
void
AverageVelocitiesMeshAdaptation<dim>::print_averaged_fields(
  const std::string &label)
{
  deallog << label << std::endl;
  deallog << "  Average velocity         : "
          << field_l2_norm(*this->average_velocities->get_average_velocities())
          << std::endl;
  deallog << "  Reynolds normal stresses : "
          << field_l2_norm(
               this->average_velocities->get_reynolds_normal_stresses())
          << std::endl;
  deallog << "  Reynolds shear stresses  : "
          << field_l2_norm(
               this->average_velocities->get_reynolds_shear_stresses())
          << std::endl;
}

template <int dim>
void
AverageVelocitiesMeshAdaptation<dim>::advance(
  const unsigned int number_of_steps)
{
  for (unsigned int step = 0; step < number_of_steps; ++step)
    {
      this->simulation_control->integrate();

      set_solution(this->simulation_control->get_current_time());

      // postprocess_fd is the production function that accumulates the
      // averaged fields, which guarantees that the averaging is driven exactly
      // as it is during a simulation.
      this->postprocess_fd(false);
      this->finish_time_step();
    }
}

template <int dim>
void
AverageVelocitiesMeshAdaptation<dim>::run()
{
  GridGenerator::hyper_cube(*this->triangulation, 0, 1);
  this->triangulation->refine_global(3);

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

  this->simulation_parameters.physical_properties_manager.initialize(
    physical_properties);

  this->setup_dofs();

  // Accumulate the averaged fields on the initial mesh.
  advance(5);
  print_averaged_fields("Averaged fields before the mesh refinement");

  // Refine the mesh through the production path. The averaged fields are
  // transferred by AverageVelocities::prepare_for_mesh_adaptation and
  // AverageVelocities::post_mesh_adaptation, which are called by
  // refine_mesh_uniform through transfer_solution.
  this->refine_mesh_uniform();
  print_averaged_fields("Averaged fields after the mesh refinement");

  // Keep accumulating on the refined mesh. If the transfer had lost the
  // accumulators, the averaged fields would restart from zero here while the
  // total averaging time would keep spanning the whole averaging window.
  advance(5);
  print_averaged_fields("Averaged fields at the end of the simulation");
}

void
test()
{
  ParameterHandler              prm;
  SimulationParameters<2>       NSparam;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 0;

  NSparam.declare(prm, size_of_subsections);
  NSparam.parse(prm);

  // Manually alter some of the default parameters of the solver
  NSparam.simulation_control.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  NSparam.simulation_control.dt       = 0.1;
  NSparam.simulation_control.time_end = 1.0;

  // A second order element is used so that the reynolds stresses, which are
  // quadratic in space, belong to the finite element space and are therefore
  // independent of the mesh. The matrix-free solver requires an equal order
  // interpolation.
  NSparam.fem_parameters.velocity_degree = 2;
  NSparam.fem_parameters.pressure_degree = 2;

  // The averaging starts with the first time step so that all the time steps
  // of the simulation contribute to the averaged fields.
  NSparam.post_processing.calculate_average_velocities        = true;
  NSparam.post_processing.initial_time_for_average_velocities = 0.0;

  // The uniform refinement is triggered manually by the test, but the maximum
  // refinement level must allow it.
  NSparam.mesh_adaptation.type = Parameters::MeshAdaptation::Type::uniform;
  NSparam.mesh_adaptation.maximum_refinement_level = 10;

  NSparam.linear_solver.at(PhysicsID::fluid_dynamics).verbosity =
    Parameters::Verbosity::quiet;
  NSparam.physics_solving_strategy.at(PhysicsID::fluid_dynamics).verbosity =
    Parameters::Verbosity::quiet;
  NSparam.boundary_conditions.createNoSlip();

  AverageVelocitiesMeshAdaptation<2> problem_2d(NSparam);
  problem_2d.run();
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
