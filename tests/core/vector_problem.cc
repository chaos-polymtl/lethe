/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

/**
 * @brief This code tests averaging values in time with Trilinos vectors.
 */

// Deal.II includes
#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/component_mask.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/postprocessing_velocities.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  Triangulation<3> tria(
    typename Triangulation<3>::MeshSmoothing(
      Triangulation<3>::smoothing_on_refinement |
      Triangulation<3>::smoothing_on_coarsening));
      
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<3> dof_handler(tria);
    
  FESystem<3> fe(FE_Q<3>(1), 3, FE_Q<3>(1), 1);
  
  dof_handler.distribute_dofs(fe);
  
  Vector<double> dummy_solution;
  
  dummy_solution.reinit(dof_handler.n_dofs());
  
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  
  Vector<double> cell_dummy_solution(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
        cell_dummy_solution = 0;
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const auto comp_i = fe.system_to_component_index(i).first;

          cell_dummy_solution(i) = 1.0*(float(comp_i)+1);
          
        }
        
       cell->get_dof_indices(local_dof_indices);
       
       for (unsigned int i = 0; i < dofs_per_cell; ++i)
          dummy_solution(local_dof_indices[i]) += cell_dummy_solution(i);
      
    }
    
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(3);
  
  ComponentMask velocity_mask = fe.component_mask(velocities);
  ComponentMask pressure_mask = fe.component_mask(pressure);
  
  
  std::vector<IndexSet> index_set_velocity = DoFTools::locally_owned_dofs_per_component(dof_handler, velocity_mask); 
  std::vector<IndexSet> index_set_pressure = DoFTools::locally_owned_dofs_per_component(dof_handler, pressure_mask); 

  // Vector<double> velocity_solution(dof_handler.n_locally_owned_dofs());	
  // Vector<double> pressure_solution(dof_handler.n_locally_owned_dofs());
  
  double correction_norm = 0.0;
  double max_correction = DBL_MIN;
  
  for (unsigned int d = 0; d < 3; ++d)
  {  
    
    for (auto j = index_set_velocity[d].begin(); j !=index_set_velocity[d].end(); j++)
    {
      correction_norm += dummy_solution[*j]*dummy_solution[*j];
      
      if (dummy_solution[*j]>max_correction)
      {
        max_correction = dummy_solution[*j];
      }
    }
  }
  
  deallog << "||u||_L2 : " << std::sqrt(correction_norm) << std::endl;
  deallog << "||u||_Linfty : " << max_correction << std::endl;
  
  correction_norm = 0.0;
  max_correction = DBL_MIN;
  
  for (auto j = index_set_pressure[3].begin(); j !=index_set_pressure[3].end(); j++)
    {
      correction_norm += dummy_solution[*j]*dummy_solution[*j];
      
      if (dummy_solution[*j]>max_correction)
      {
        max_correction = dummy_solution[*j];
      }
    }
  
  deallog << "||p||_L2 : " << std::sqrt(correction_norm) << std::endl;
  deallog << "||p||_Linfty : " << max_correction << std::endl;

}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
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
