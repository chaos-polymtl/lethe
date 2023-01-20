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
 * This class compartmentalizes the geometry/reactor based on physical
 * values such as velocity, magnetic field, etc., and calculates the
 * flux between each compartment based on the velocity.
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
   * @brief Set up grid and compartmentalize first based on electric values
   * and then find the final pattern based on velocity
   */
  void
  test();

private:
  /**
   * @brief Read the physical property (in this problem, the electric field)
   * and write it into a vector of a vector.
   */
  auto
  read_electric_field() -> std::vector<std::vector<double>>;

  /**
   * @brief Read the velocity magnitude and write it into a map.
   */
  auto
  read_velocity_magnitude()
    -> std::map<typename Triangulation<dim>::active_cell_iterator, double>;

  /**
   * @brief Generate subdivided_cylinder mesh
   */
  void
  generate_cylindrical_grid();

  /**
   * @brief 1) Sorts the cell based on their
   *  values of the physical properties
   *  (average physical values, in particular, electric field)
   *  2) agglomerates cell within tolerance and then 3) break the
   * clusters that are not connected but are in one group
   * because of their tolerance.
   */
  auto
  sort_agglomeration_deagglomeration_emw()
    -> std::map<int,
                std::vector<typename Triangulation<dim>::active_cell_iterator>>;

  /**
   * @brief this function gets the first set of compartments based on one
   * of the physical values (for example, electric field) from
   * "sort_agglomeration_deagglomeration_emw()" for each compartment
   * this function treats each compartment as a new problem to be
   * compartmentalized and based on the second physical property (velocity) it
   * 1) sorts the cell, 2) agglomerate within the tolerance, 3) deagglomerate
   * the not connected parts
   */
  void
  overlaid_map();

  /**
   * @brief Write the pvd file of the set of compartments based on first
   * physical property
   */
  void
  write_file_compartments_first_field(Vector<double> &    compartments_final,
                                      const double &      time,
                                      const unsigned int &step_number);

  /**
   * @brief Write the pvd file of the final set of compartments
   */
  void
  write_ultimate_file(Vector<double> &    compartments_ultimate,
                      const double &      time,
                      const unsigned int &step_number);


  // Vector of vector to store cell_index and the average physical value of the
  // cell
  std::vector<std::vector<double>> primer_matrix_index_emw;
  // Vector of vector to store cell_index and the average velocity magnitude of
  // the cell
  std::map<typename Triangulation<dim>::active_cell_iterator, double>
    primer_matrix_index_velocity;
  // Map of final set to store the compartment_id and the cells that each
  // compartment includes
  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    final_set;

  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    overlaid_set;

  std::map<int, std::vector<typename Triangulation<dim>::active_cell_iterator>>
    ultimate_map;


  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::vector<double>>
    map_of_all_cells_and_velocity_vectors;
  // Map of paris of cell and cell_index
  std::map<int, typename Triangulation<dim>::active_cell_iterator>
                                            cells_and_indices;
  TimerOutput                               computing_timer;
  PVDHandler                                grid_pvdhandler;
  MPI_Comm                                  mpi_communicator;
  parallel::distributed::Triangulation<dim> triangulation;
  CPCalculatingParameters                   cp_parameters;
};



#endif // lethe_compartmentalization_h
