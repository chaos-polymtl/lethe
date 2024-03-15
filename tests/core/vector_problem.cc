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
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  
  Triangulation<3> tria(
    typename Triangulation<3>::MeshSmoothing(
      Triangulation<3>::smoothing_on_refinement |
      Triangulation<3>::smoothing_on_coarsening));
      
  GridGenerator::hyper_cube(tria, -1, 1);
  DoFHandler<3> dof_handler(tria);
    
  FESystem<3> fe(FE_Q<3>(2), 3, FE_Q<3>(1), 1);
  
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
          
          if (comp_i < 3)
          {
            cell_dummy_solution(i) = 1.0*(float(comp_i)+1);
          }
        }
        
       cell->get_dof_indices(local_dof_indices);
       
       for (unsigned int i = 0; i < dofs_per_cell; ++i)
          dummy_solution(local_dof_indices[i]) += cell_dummy_solution(i);
      
    }
    
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(3);
  
  ComponentMask velocity_mask = fe.component_mask(velocities);
  ComponentMask pressure_mask = fe.component_mask(pressure);
  
  
  std::vector<IndexSet> index_set = DoFTools::locally_owned_dofs_per_component(dof_handler, velocity_mask); 
  
  Vector<double> velocity_solution(3*index_set[0].n_elements());	
  
  unsigned int k = 0;
  for (unsigned int d = 0; d < 3; ++d)
  {
    
    unsigned int index_set_size = index_set[d].n_elements();
    
    std::cout << "IndexSet comp " << d << " size = " << index_set_size << std::endl;
    
    for (auto j = index_set[d].begin(); j !=index_set[d].end(); j++, k++)
    {
      velocity_solution[k] = dummy_solution[*j];
      
      std::cout << "DoF index = " <<  *j << " dummy_velocity = " << velocity_solution[k] <<'\n';
    }
  }
  
  for (unsigned int i = 0; i < 3*index_set[0].n_elements(); ++i)
  {
    auto comp_i = i/index_set[0].n_elements();
    
    std::cout << "velocity comp_i " <<  comp_i << " = " << velocity_solution[i] <<'\n';
  }
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
