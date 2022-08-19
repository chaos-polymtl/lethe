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
  RPTCellReconstruction(
    Parameters::RPTParameters &              rpt_parameters,
    Parameters::RPTReconstructionParameters &rpt_reconstruction_parameters,
    Parameters::DetectorParameters &         rpt_detector_parameters);

  /**
   * @brief Set up grid and find positions for every unknown particle positions.
   */
  void
  execute_cell_reconstruction();

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
   *
   * @param id Index of the unknown particle, useful when "positions_evaluation" is enable
   */
  void
  find_unknown_position(std::vector<double> &particle_reconstruction_counts,
                        unsigned int         id);

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
   * @param particle_reconstruction_counts Counts for every detectors associated
   *                                       to the unknown position of the
   * particle
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
   * @brief Export positions and cell volumes in .csv or .dat
   */
  void
  export_positions();

  /**
   * @brief Generate files for visualization (pvtu/vtu)
   */
  void
  visualize_positions();


private:
  Triangulation<dim>                      triangulation;
  Parameters::RPTParameters               parameters;
  Parameters::RPTReconstructionParameters reconstruction_parameters;
  Parameters::DetectorParameters          detector_parameters;
  std::vector<Detector<dim>>              detectors;

  // All counts of the unknown particle positions
  // Data structure : [[detector_counts_0_particle_0, ...,
  //                    detector_counts_n_particle_0],
  //                   [...],
  //                   [detector_counts_0_particle_n, ...,
  //                   detector_counts_n_particle_n]]
  std::vector<std::vector<double>> reconstruction_counts;

  // Map of all calculated counts of vertices (key : vertex_id)
  // Data structure : {{vertex_id_0, <vertex_position_0(x, y, z),
  //                    [detector_counts_0_particle_0, ...]>},
  //                   {...},
  //                   {vertex_id_m, <vertex_position_m(x, y, z),
  //                   [detector_counts_0_particle_0, ...]>}}
  std::map<unsigned int, std::pair<Point<dim>, std::vector<double>>>
    map_vertices_index;

  std::vector<Point<dim>> reconstruction_positions; // Found positions
  std::vector<double>
                            cells_volumes; // Cell volumes of the cell that contains positions
  std::vector<Point<dim>>   known_positions;  // Known positions to analyse
  std::vector<unsigned int> final_cell_level; // Level of the best cells

  // Status of the best cells after analyse of reconstructed positions vs real
  // positions (many status can be associated to a cell)
  //  - right_cell :    The reconstructed position is in the same cell than the
  //  real
  //                    position
  //  - cost_function : There was many cell candidates at refined mesh and a
  //                    cost function was needed to get the best cell
  //  - parent_cell :   The best cell is not at the highest refined level
  //  - wrong_cell :    The reconstructed position is not in the same cell than
  //                    the real position
  std::vector<std::string> cell_status;

  TimerOutput computing_timer;
};


#endif // lethe_rpt_cell_reconstruction_h
