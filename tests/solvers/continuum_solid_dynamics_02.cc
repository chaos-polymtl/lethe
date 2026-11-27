// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the initial and boundary conditions of the continuum
 * solid dynamics. A constant initial condition is set on the unit hyper cube
 * and the solution is evaluated:
 * - inside the domain, where it equals the initial condition;
 * - on a noslip boundary, where the solid velocity is zero and the solid
 *   volume fraction is free;
 * - on a function boundary, where the solid velocity and the solid volume
 *   fraction are imposed;
 * - on a slip boundary, where only the normal solid velocity is zero.
 */

#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

#include <../tests/tests.h>

#include <memory>
#include <string>

using namespace dealii;

/**
 * @brief Print the solid velocity and solid volume fraction at a point.
 *
 * @param[in] label Description of the point.
 * @param[in] multiphysics Multiphysics interface holding the solution.
 * @param[in] point Point at which the solution is evaluated.
 */
template <int dim>
void
print_value_at_point(const std::string          &label,
                     MultiphysicsInterface<dim> &multiphysics,
                     const Point<dim>           &point)
{
  Vector<double> values(dim + 1);
  VectorTools::point_value(
    multiphysics.get_mapping(PhysicsID::continuum_solid_dynamics),
    multiphysics.get_dof_handler(PhysicsID::continuum_solid_dynamics),
    multiphysics.get_solution(PhysicsID::continuum_solid_dynamics),
    point,
    values);

  deallog << label << ":";
  for (const double value : values)
    deallog << " " << value;
}

template <int dim>
void
test()
{
  // Boundary ids of the colorized hyper cube: 0 (x = 0), 1 (x = 1), 2 (y = 0)
  auto tria =
    std::make_shared<parallel::distributed::Triangulation<dim>>(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(*tria, 0., 1., true);
  tria->refine_global(2);

  SimulationParameters<dim>     solver_parameters;
  ParameterHandler              prm;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 3;
  size_of_subsections.manifolds           = 1;
  solver_parameters.declare(prm, size_of_subsections);

  // Initial condition: u_s = 1, v_s = 2, (w_s = 3), alpha = 0.4
  const std::string initial_condition =
    (dim == 2) ? "1; 2; 0.4" : "1; 2; 3; 0.4";
  // Function boundary: u_s = 0.5, other velocity components 0, alpha = 0.3
  const std::string boundary_function =
    (dim == 2) ? "0.5; 0; 0.3" : "0.5; 0; 0; 0.3";

  prm.enter_subsection("initial conditions");
  prm.enter_subsection("continuum solid dynamics");
  prm.set("Function expression", initial_condition);
  prm.leave_subsection();
  prm.leave_subsection();

  prm.enter_subsection("boundary conditions continuum solid dynamics");
  prm.set("number", "3");
  prm.enter_subsection("bc 0");
  prm.set("id", "0");
  prm.set("type", "noslip");
  prm.leave_subsection();
  prm.enter_subsection("bc 1");
  prm.set("id", "1");
  prm.set("type", "function");
  prm.enter_subsection("function");
  prm.set("Function expression", boundary_function);
  prm.leave_subsection();
  prm.leave_subsection();
  prm.enter_subsection("bc 2");
  prm.set("id", "2");
  prm.set("type", "slip");
  prm.leave_subsection();
  prm.leave_subsection();

  solver_parameters.parse(prm);
  solver_parameters.multiphysics.fluid_dynamics           = true;
  solver_parameters.multiphysics.continuum_solid_dynamics = true;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                          tria,
                                          simulation_control,
                                          pcout,
                                          FluidDynamicsFormulation::VANS);
  multiphysics.setup_dofs();
  multiphysics.set_initial_conditions();

  // Mid-points of the domain and of the boundaries (away from the corners)
  Point<dim> center, noslip_point, function_point, slip_point;
  for (unsigned int d = 0; d < Point<dim>::dimension; ++d)
    {
      center[d]         = 0.5;
      noslip_point[d]   = 0.5;
      function_point[d] = 0.5;
      slip_point[d]     = 0.5;
    }
  noslip_point[0]   = 0.;
  function_point[0] = 1.;
  slip_point[1]     = 0.;

  print_value_at_point("Initial condition (center)", multiphysics, center);
  print_value_at_point("Noslip boundary", multiphysics, noslip_point);
  print_value_at_point("Function boundary", multiphysics, function_point);
  print_value_at_point("Slip boundary", multiphysics, slip_point);
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
      test<3>();
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
