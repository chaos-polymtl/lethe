/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
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
 */

#ifndef lethe_compartmentalization_h
#define lethe_compartmentalization_h

/**
 * Compartmentalizes a reactor geometry based on physical values such as
 * velocity, magnetic field, etc. It calculates the flux between each
 * compartment based on the velocity between compartments.
 */

#include <deal.II/base/config.h>

#include <core/pvd_handler.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <rpt/parameters_cp.h>

#include <map>
#include <unordered_set>
#include <vector>

using namespace dealii;

template <int dim>
class Compartmentalization
{
public:
  /**
   * @brief Constructor for the Compartmentalization
   */
  Compartmentalization<dim>(CPCalculatingParameters &CPparameters);

  /**
   * @brief Compartmentalizes the cells based on electric field values
   * and determines the overlaid pattern using cell velocity magnitude
   * print out the volume of each compartment
   * print out the average of electric field at each compartment
   * print out the matrix of inlet/outlet flux for the cylinder and flux between
   * compartments
   */
  void
  test();

private:
  /**
   * @brief Read the electric field at each cell, derived from COMSOL,
   * and write the pair of cell_index and electric field into a vector of a
   * vector to use it for compartmentalization
   */
  std::vector<std::vector<double>>
  read_electric_field_in_vector();

  /**
   * @brief Read the electric field at each cell, derived from COMSOL
   * and write each cell as key and the associated electric field in a map
   * to use it for calculating the average electric field at the final
   * compartments
   */
  std::map<typename Triangulation<dim>::active_cell_iterator, double>
  read_electric_field_in_map();


  /**
   * @brief Read the velocity magnitude and write the pair of cell and associated velocity to that cell into a map
   */
  std::map<typename Triangulation<dim>::active_cell_iterator, double>
  read_velocity_magnitude();

  /**
   * @brief Read the velocity vector at the center point of cells and write it into a vector
   * add the cell as the key and associated velocity vector to a map
   */
  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::vector<double>>
  read_velocity_vector();


  /**
   * @brief Generate subdivided_cylinder mesh
   */
  void
  generate_cylindrical_grid();
  /**
   * @brief output the cell center for pyvista
   */
  void
  write_cell_center();

  /**
   * @brief First, sorts the cell based on their values of the electric field.
   * Then, agglomerates cells in groups according to tolerances and
   * breaks the clusters accordingly.
   */
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
  compartmentalize_first_step();

  /**
   * @brief This function performs a two-step compartmentalization process. Firstly,
   * it obtains the initial set of compartments by analyzing the electric field
   * values of each cell, utilizing the
   * compartmentalize_first_step()" function. Subsequently, each
   * compartment is treated as a separate problem to be compartmentalized
   * further. The second physical property, velocity, is employed to carry out
   * the following steps within each compartment: 1) sorts the cell, 2)
   * agglomerate within the tolerance, 3) deagglomerate the not connected parts
   */
  void
  overlaid_map();

  /**
   * @brief Write the pvd file of the set of compartments based on electric field
   */
  void
  write_file_compartments_first_field(Vector<double> &    compartments_final,
                                      const double &      time,
                                      const unsigned int &step_number);

  /**
   * @brief Write the pvd file of the final set of compartments (overlaid map)
   */
  void
  write_ultimate_file(Vector<double> &    compartments_ultimate,
                      const double &      time,
                      const unsigned int &step_number);

  /**
   * @brief Calculate the flux between the compartments based on their shared
   * surface and velocity
   */
  void
  flux_calculation_between_compartments();

  /**
   * @brief Calculate the inlet at outlet flux at the two boundaries
   */
  void
  inlet_outlet_flux();

  /**
   * @brief Prepare a matrix of fluxes (inlet,outlet,interconnecting fluxes) inorder to use as the input to python code
   * write it in a txt file
   */
  void
  output_flux();

  // Map contains cell as the key and the associated value of the electric field
  // at each cell
  std::map<typename Triangulation<dim>::active_cell_iterator, double>
    map_cell_electric_field;

  // Vector of vector to store cell_index and the average value of electric
  // field in each cell
  std::vector<std::vector<double>> matrix_index_electric_field;

  // Vector of vector to store cell_index and the average velocity magnitude of
  // the cell
  std::map<typename Triangulation<dim>::active_cell_iterator, double>
    matrix_index_velocity;

  // Map of cell (active_cell_iterator) as key and the associated velocity
  // vector
  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::vector<double>>
    map_of_all_cells_and_velocity_vectors;

  // Map of final set to store the compartment_id and the cells that each
  // compartment includes
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    final_set;

  // Define the map that will contain the elements for flux calculation
  std::map<std::string, double> flux;

  // Map containing the final compartment_ids as key and a vector of cell
  // inside each compartment
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    overlaid_set;

  // Vector of inlet flux
  std::map<std::string, double> inlet_flux;

  // Vector of outlet flux
  std::map<std::string, double> outlet_flux;

  // Map of paris of cell_index (int) and cell (active_cell_iterator)
  std::map<int, typename Triangulation<dim>::active_cell_iterator>
    cells_and_indices;

  TimerOutput                               computing_timer;
  PVDHandler                                grid_pvdhandler;
  MPI_Comm                                  mpi_communicator;
  parallel::distributed::Triangulation<dim> triangulation;
  CPCalculatingParameters                   cp_parameters;
};



#endif // lethe_compartmentalization_h
