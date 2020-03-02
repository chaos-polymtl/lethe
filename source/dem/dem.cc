/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Shahab Golshan, Polytechnique Montreal, 2019-
 */

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <dem/dem.h>
#include <dem/dem_properties.h>

#include <fstream>
#include <iostream>

template <int dim>
DEMSolver<dim>::DEMSolver(DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
  , parameters(dem_parameters)
  , triangulation(this->mpi_communicator)
  , property_pool(DEM::get_number_properties())
  , mapping(1)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , particle_handler(triangulation, mapping, DEM::get_number_properties())
{}

template <int dim>
void
DEMSolver<dim>::read_mesh()
{
  // GMSH input
  if (parameters.mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);
      std::ifstream input_file(parameters.mesh.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (parameters.mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        triangulation,
        parameters.mesh.grid_type,
        parameters.mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  const int initialSize = parameters.mesh.initialRefinement;
  triangulation.refine_global(initialSize);
}

template <int dim>
void
DEMSolver<dim>::solve()
{
  // Reading mesh
  read_mesh();

  // Initializing variables
  int    DEM_step        = 0;
  double DEM_time        = 0;
  int    number_of_steps = parameters.simulationControl.final_time_step;
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
                                                          cell_neighbor_list;
  std::vector<std::map<int, pp_contact_info_struct<dim>>> pairs_in_contact_info(
    parameters.simulationControl.total_particle_number);
  std::vector<std::map<int, pw_contact_info_struct<dim>>> pw_pairs_in_contact(
    parameters.simulationControl.total_particle_number);
  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;
  DEM::DEMProperties<dim>                      properties_class;
  Tensor<1, dim>                               g;
  if (dim == 3)
    {
      g[0] = parameters.physicalProperties.gx;
      g[1] = parameters.physicalProperties.gy;
      g[2] = parameters.physicalProperties.gz;
    }
  if (dim == 2)
    {
      g[0] = parameters.physicalProperties.gx;
      g[1] = parameters.physicalProperties.gy;
    }

  // Initilization of classes and building objects
  DEM_iterator<dim>             iterator_object;
  VelocityVerletIntegrator<dim> integrator_object;
  // ***** I need to choose the contact model based on input file
  PPBroadSearch<dim> pp_broad_search_object;
  PPFineSearch<dim>  pp_fine_search_object;
  // PPLinearForce<dim> pp_force_object;
  PPNonLinearForce<dim> pp_force_object;
  // PWLinearForce<dim> pw_force_object;
  PWNonLinearForce<dim> pw_force_object;
  PWBroadSearch<dim>    pw_broad_search_object;
  PWFineSearch<dim>     pw_fine_search_object;

  std::vector<std::pair<std::string, int>> properties =
    properties_class.get_properties_name();

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbors_object;
  cell_neighbor_list = cell_neighbors_object.find_cell_neighbors(triangulation);

  // Finding boundary cells
  FindBoundaryCellsInformation<dim> boundary_cell_object;
  boundary_cells_information =
    boundary_cell_object.find_boundary_cells_information(triangulation);

  // DEM engine iterator:
  while (DEM_step < number_of_steps)
    {
      iterator_object.engine(particle_handler,
                             triangulation,
                             DEM_step,
                             DEM_time,
                             cell_neighbor_list,
                             pairs_in_contact_info,
                             boundary_cells_information,
                             pw_pairs_in_contact,
                             parameters,
                             g,
                             properties,
                             property_pool,
                             &pp_force_object,
                             &pw_force_object,
                             &integrator_object,
                             &pp_broad_search_object,
                             &pp_fine_search_object,
                             &pw_broad_search_object,
                             &pw_fine_search_object,
                             computing_timer);
    }

  while (parameters.simulation_control.integrate())
    {
      printTime(this->pcout, parameters.simulation_control);
    }
}

template class DEMSolver<2>;
template class DEMSolver<3>;
