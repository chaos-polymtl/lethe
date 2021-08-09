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

*
* Authors: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

#ifndef lethe_rpt_cell_reconstruction_h
#define lethe_rpt_cell_reconstruction_h

/**
 * This class allows cell position reconstruction of unknown particle positions
 * with counts value for many detector. It finds the center position of the best
 * cell that its vertices counts range contains the counts value of unknown
 * particle positions.
 */

#include <deal.II/base/config.h>

#include <deal.II/base/timer.h>

#include <deal.II/grid/tria.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_visualization.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

#include <map>
#include <vector>

using namespace dealii;

template <int dim>
class RPTCellReconstruction
{
public:
  /**
   * @brief Constructor for the RPTCellReconstruction
   *
   * @param detectors Vector of all detector object
   *
   * @param rpt_parameters All parameters for the RPT
   *
   */
  RPTCellReconstruction(RPTCalculatingParameters &rpt_parameters);

  /**
   * @brief Set up grid and find positions for every unknown particle positions.
   */
  void
  execute_cell_reconstruction();

  /**
   * @brief Create grid of the reactor vessel (cylinder).
   */
  void
  create_grid();

  /**
   * @brief Read text file for detector positions and assign them to detector
   * objects
   */
  void
  assign_detector_positions();

  /**
   * @brief Calculate counts at vertices for the coarse mesh.
   */
  void
  set_coarse_mesh_counts();

  /**
   * @brief Find the position of the best cell that contains the unknown particle
   * position
   *
   * @param particle_reconstruction_counts Counts the the current unknown particle
   * position for every detector
   */
  void
  find_unknown_position(std::vector<double> &particle_reconstruction_counts);

  void
  find_all_positions();

  /**
   * @brief Find the vertices position of coarse mesh vertices or the candidate
   * cells
   *
   * @param level Level of the current cells of the mesh
   *
   * @param parent_cell_indexes List of the parent cells which are candidate.
   *                            Is {-1} by default when executing coarse mesh
   *                            (no parent cells)
   */
  void
  find_vertices_positions(unsigned int     level,
                          std::vector<int> parent_cell_indexes = {-1});
  /**
   * @brief Find the cell candidates indexes for the unknown particle positions
   * cells
   *
   * @param level Level of the current cells of the mesh
   *
   * @param parent_cell_indexes List of the parent cells which are candidate.
   *                            Is {-1} by default when executing coarse mesh
   *                            (no parent cells)
   */
  std::vector<int>
  find_cells(unsigned int         level,
             std::vector<double> &particle_reconstruction_counts,
             std::vector<int>     parent_cell_indexes = {-1});

  /**
   * @brief Find the best cell candidates indexes for the unknown particle positions
   * with the least squared method comparing vertices counts for every detectors
   *
   * @param level Level of the current cells of the mesh
   *
   * @param candidate List of cell candidates
   */
  int
  find_best_cell(unsigned int         level,
                 std::vector<double> &particle_reconstruction_counts,
                 std::vector<int>     candidate);

  /**
   * @brief Calculate counts at new vertices
   */
  void
  calculate_counts();

  /**
   * @brief Read the .reconstruction file to extract counts of unknown particle
   * positions
   */
  void
  read_counts();

  /**
   * @brief Export positions and cell volumes in .csv or .dat
   */
  void
  export_positions();

  /**
   * @brief
   */
  void
  visualize_positions();

private:
  Triangulation<dim>         triangulation;
  RPTCalculatingParameters   parameters;
  std::vector<Detector<dim>> detectors;

  std::vector<std::vector<double>>
    reconstruction_counts; // All counts of the unknown particle positions
  std::map<unsigned int, std::pair<Point<dim>, std::vector<double>>>
    map_vertices_index; // Map of all calculated counts of vertices
  std::vector<Point<dim>> reconstruction_positions; // Positions found
  std::vector<double>
    cells_volumes; // Cell volumes of the cell that contains positions

  TimerOutput computing_timer;
};


#endif // lethe_rpt_cell_reconstruction_h
