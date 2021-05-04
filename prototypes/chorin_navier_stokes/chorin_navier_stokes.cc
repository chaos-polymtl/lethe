/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

//#include "boundaryconditions.h"
//#include "exactsolutions.h"
//#include "forcingfunctions.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <fstream>
#include <cmath>

// Finally, this is as in previous programs:
using namespace dealii;

enum SimulationCases
{
  MMS           = 0,
  TaylorCouette = 1,
};

template <int dim>
class ChorinNavierStokes
{
public:
  ChorinNavierStokes(const unsigned int degreeVelocity,
                     const unsigned int degreePressure);
  ~ChorinNavierStokes();
  void
  runMMS();
  void
  runCouette();
  
  void
  run();

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

private:
  void
  make_cube_grid(int refinementLevel);
  void
  refine_grid();
  void
  refine_mesh();
  void
  refine_mesh_uniform();
  void
  setup_dofs();
  void
  initialize_system();

  void
  calculateL2Error();
  void
  output_results(const unsigned int cycle) const;

  double             viscosity;
  Triangulation<dim> triangulation;

  // FE system for velocity
  FESystem<dim>      fe_velocity;
  DoFHandler<dim>    dof_handler_velocity;
  SparsityPattern    sparsity_pattern_velocity;

  // FE system for pressure
  FESystem<dim>      fe_pressure;
  DoFHandler<dim>    dof_handler_pressure;
  SparsityPattern    sparsity_pattern_pressure;

};


// Constructor
template <int dim>
ChorinNavierStokes<dim>::ChorinNavierStokes(const unsigned int degreeVelocity,
                                            const unsigned int degreePressure)
  : viscosity(1)

  // Initialise FE system for velocity
  , fe_velocity(FE_Q<dim>(degreeVelocity), dim)
  , dof_handler_velocity(triangulation)

  // Initialise FE system for pressure
  , fe_pressure(FE_Q<dim>(degreePressure), 1)
  , dof_handler_pressure(triangulation)
{}

template <int dim>
ChorinNavierStokes<dim>::~ChorinNavierStokes()
{
  triangulation.clear();
}

template <int dim>
void
ChorinNavierStokes<dim>::make_cube_grid(int refinementLevel)
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(refinementLevel);
}

template <int dim>
void
ChorinNavierStokes<dim>::refine_grid()
{
  triangulation.refine_global(1);
}

template <int dim>
void
ChorinNavierStokes<dim>::setup_dofs()
{
  // Distribute DOFs
  dof_handler_velocity.distribute_dofs(fe_velocity);
  dof_handler_pressure.distribute_dofs(fe_pressure);

  // Output information
  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Number of velocity degrees of freedom: " << dof_handler_velocity.n_dofs()
            << std::endl
            << "   Number of pressure degrees of freedom: " << dof_handler_pressure.n_dofs()
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::initialize_system()
{
  // Build sparsity patterns
  DynamicSparsityPattern dsp_velocity(dof_handler_velocity.n_dofs(), dof_handler_velocity.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp_velocity);
  sparsity_pattern_velocity.copy_from(dsp_velocity);

  DynamicSparsityPattern dsp_pressure(dof_handler_pressure.n_dofs(), dof_handler_pressure.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp_pressure);
  sparsity_pattern_pressure.copy_from(dsp_pressure);

  std::ofstream sparsity_velocity_out("sparsity_pattern_velocity.svg");
  sparsity_pattern_velocity.print_svg(sparsity_velocity_out);

  std::ofstream sparsity_pressure_out("sparsity_pattern_pressure.svg");
  sparsity_pattern_pressure.print_svg(sparsity_pressure_out);
}

template <int dim>
void
ChorinNavierStokes<dim>::refine_mesh()
{}

template <int dim>
void
ChorinNavierStokes<dim>::refine_mesh_uniform()
{}

template <int dim>
void
ChorinNavierStokes<dim>::output_results(const unsigned int cycle) const
{
  /*
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler_pressure);
  //data_out.add_data_vector(solution, "pressure");

  data_out.build_patches();

  std::ofstream output("output_pressure-" + Utilities::int_to_string(cycle, 4) + ".vtu");
  data_out.write_vtu(output);
  */
  
  std::ofstream out("grid.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid.svg" << std::endl;
  
}

// Find the l2 norm of the error between the finite element solution and the exact solution
template <int dim>
void
ChorinNavierStokes<dim>::calculateL2Error()
{}

// Run function while building to test
template <int dim>
void
ChorinNavierStokes<dim>::run()
{
  make_cube_grid(3);
  refine_grid();
  setup_dofs();
  initialize_system();
  output_results(0);
}

template <int dim>
void
ChorinNavierStokes<dim>::runMMS()
{}


template <int dim>
void
ChorinNavierStokes<dim>::runCouette()
{}

int
main()
{
  try
    {
      ChorinNavierStokes<2> problem_2d(1, 1); // degreeVelocity, degreePressure
      problem_2d.run();
      //problem_2d.runCouette();
      //problem_2d.runMMS();
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