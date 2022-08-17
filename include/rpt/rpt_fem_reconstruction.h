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

/**
 * This class calculates an L2 projection of the counts on a given
 * triangulation then uses this projection to reconstruct
 * the position of the particles
 */

#ifndef lethe_rpt_fem_reconstruction_h
#define lethe_rpt_fem_reconstruction_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>
#include <rpt/rpt_utilities.h>


using namespace dealii;

template <int dim>
class RPTFEMReconstruction
{
public:
  /**
   * @brief Constructor for the RPTFEMReconstruction
   *
   * @param RPTparameters All parameters and positions needed for the count
   * calculation and the position reconstruction
   */
  RPTFEMReconstruction(Parameters::RPTParameters &rpt_parameters,
                       Parameters::RPTFEMReconstructionParameters
                         &rpt_fem_reconstruction_parameters,
                       Parameters::DetectorParameters &rpt_detector_parameters)
    : fe(1)
    , dof_handler(triangulation)
    , mapping(FE_SimplexP<dim>(1))
    , parameters(rpt_parameters)
    , fem_reconstruction_parameters(rpt_fem_reconstruction_parameters)
    , detector_parameters(rpt_detector_parameters)
    , computing_timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
  {}

  /**
   * @brief Using the L2 projection, builds a dictionary of photon counts
   * recorded by a detector for given positions inside a  reactor.
   */
  void
  L2_project();

  /**
   * @brief Reconstruction of the particle positions with the counts from the
   * detectors.
   */
  void
  rpt_fem_reconstruct();


private:
  /**
   * @brief Set up the triangulation that covers the domain of the reactor.
   */
  void
  setup_triangulation();

  /**
   * @brief Set up the system to be solved by enumerating the degrees of
   * freedom and set up the matrix and vector objects to hold the system data.
   */
  void
  setup_system();

  /**
   * @brief Assemble the matrix and the right hand side that form the linear system
   * that needs to be solved for a given detector.
   *
   * @param detector_no detector_no defines the detector for which the
   * linear system is being assembled.
   */
  void
  assemble_system(unsigned detector_no);

  /**
   * @brief Solve the linear system for a given detector to get the nodal
   * counts.
   *
   * @param detector_no detector_no defines the detector for which the linear
   * system is being solved.
   */
  void
  solve_linear_system(unsigned detector_no);

  /**
   * @brief Outputs the solution of the linear system (nodal counts) in a
   * ".vtu" file format.
   */
  void
  output_results();

  /**
   * @brief Outputs the position of each vertex of every cell of the grid and
   * the photon count at that position for every detector in a ".dat" file
   * format.
   */
  void
  output_raw_results();

  /**
   * @brief Calculates the cost with the selected cost function of the found
   * position and returns it.
   *
   * @param cell Iterator on active cells of the grid.
   *
   * @param reference_location Location of the particle in reference
   * coordinates
   *
   * @param last_constraint Constraint on the particle's location
   *
   * @param experimental_count Experimental counts from all detectors for a
   * given particle position
   *
   * @return Cost It is calculated with the cost function mentioned in the
   * parameter file.
   */
  double
  calculate_cost(
    const TriaActiveIterator<DoFCellAccessor<dim, dim, false>> &cell,
    Vector<double> &     reference_location,
    const double &       last_constraint,
    std::vector<double> &experimental_count);

  /**
   * @brief Evaluates the error vector of the found location in reference
   * coordinates and returns the error vector's L2 norm
   *
   * @param reference_location Location of the point in reference coordinates
   *
   * @param last_constraint 4th constraint on the location of the particle in
   * reference coordinates
   *
   * @return norm_error_coordinates
   */
  double
  calculate_reference_location_error(Vector<double> &reference_location,
                                     const double &  last_constraint);

  /**
   * @brief Finds the position of the particle.
   *
   * @param experimental_count experimental_count contains the experimental
   * counts of every detector for a given position.
   *
   * @param tol_reference_location tolerance on the validity of the found
   * location in the reference space
   *
   * @return 'true' if the particle's position was found and 'false' if the
   * particle's position couldn't be found
   */
  bool
  find_position_global_search(std::vector<double> &experimental_count,
                              const double         tol_reference_location);

  /**
   * @brief Finds the position of the particle by doing a local search around
   * the previously found position's cell's neighbors
   *
   * @param experimental_count experimental_count contains the experimental
   * counts of every detector for a given position.
   *
   * @param tol_reference_location tolerance on the validity of the found
   * location in the reference space
   *
   * @param cell cell in which the previous position was
   * found
   *
   * @return 'true' if the particle's position was found and 'false' if the
   * particle's position couldn't be found
   */
  bool
  find_position_local_search(
    std::vector<double> &experimental_count,
    const double         tol_reference_location,
    const typename DoFHandler<dim>::active_cell_iterator &cell);

  /**
   * @brief Reads the file with the experimental counts and finds the
   * particles' positions by calling the "find_cell" function.
   */
  void
  trajectory();

  /**
   * @brief Saves the built dictionary so it may be used later for the 3D particle position
   * reconstruction.
   */
  void
  checkpoint();

  /**
   * @brief Loads the previously built dictionary for the 3D FEM
   * reconstruction.
   */
  void
  load_from_checkpoint();

  /**
   * @brief Exports found particle positions in a the specified file format.
   * If not specified, it will be exported in a ".csv" file format.
   */
  void
  export_found_positions();

  Triangulation<dim>                         triangulation;
  FE_SimplexP<dim>                           fe;
  DoFHandler<dim>                            dof_handler;
  RPTCalculatingParameters                   rpt_parameters;
  MappingFE<dim>                             mapping;
  Parameters::RPTParameters                  parameters;
  Parameters::RPTFEMReconstructionParameters fem_reconstruction_parameters;
  Parameters::DetectorParameters             detector_parameters;
  unsigned int                               n_detector;

  AffineConstraints<double>                      constraints;
  SparseMatrix<double>                           system_matrix;
  SparsityPattern                                sparsity_pattern;
  Vector<double>                                 system_rhs;
  std::vector<Vector<double>>                    nodal_counts;
  std::vector<Detector<dim>>                     detectors;
  std::vector<Point<dim>>                        found_positions;
  typename DoFHandler<dim>::active_cell_iterator previous_position_cell;

  TimerOutput computing_timer;

  /**
   * @brief Contains the FEValues objects
   */
  struct AssemblyScratchData
  {
    AssemblyScratchData(const FiniteElement<dim> &fe,
                        const unsigned int        no_detector);
    AssemblyScratchData(const AssemblyScratchData &scratch_data);

    FEValues<dim> fe_values;
    unsigned int  detector_id;
  };

  /**
   * @brief Contains elements of the local linear system
   */
  struct AssemblyCopyData
  {
    FullMatrix<double>                   cell_matrix;
    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };

  /**
   * @brief Assembles local linear system to get the nodal values of a given
   * cell for a given detector
   *
   * @param cell iterator on active cells of the grid
   *
   * @param scratch_data contains FEValues objects
   *
   * @param copy_data contains elements of the local linear system
   */
  void
  assemble_local_system(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    AssemblyScratchData &                                 scratch_data,
    AssemblyCopyData &                                    copy_data);

  /**
   * @brief Copies local linear system elements to the global linear system
   *
   * @param copy_data contains elements of the local linear system
   */
  void
  copy_local_to_global(const AssemblyCopyData &copy_data);
};

/**
 * @brief Sets-up, assembles and solves the linear system to find the location
 * of a particle in reference coordinates.
 *
 * @param vertex_count contains the photon counts of all detectors for a s
 * elected vertex of cell.
 *
 * @param experimental_count contains the experimental counts from all
 * detectors for a given particle position.
 *
 * @param cost_function_type cost function type used to find the position of
 * the particle in the reactor
 *
 * @return Calculated reference position
 */
template <int dim>
Vector<double>
assemble_matrix_and_rhs(
  std::vector<std::vector<double>> &vertex_count,
  std::vector<double> &             experimental_count,
  Parameters::RPTFEMReconstructionParameters::FEMCostFunction
    &cost_function_type);

#endif
