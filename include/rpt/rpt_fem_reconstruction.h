/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2022 - by the Lethe authors
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
 */

/**
 * This class calculates an L2 projection of the counts on a given
 * triangulation then uses this projection to reconstruct
 * the position of the particles
 */

#ifndef lethe_rpt_fem_reconstruction_h
#define lethe_rpt_fem_reconstruction_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/particles/particle_handler.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/particle_detector_interactions.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>
#include <rpt/rpt_utilities.h>


#define FORCE_USE_OF_TRILINOS
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

using namespace dealii;


template <int dim>
class RPTL2Projection
{
public:
  /**
   * @brief Constructor for the RPTL2Projection
   *
   * @param rpt_parameters General parameters for photon count calculation with the Monte Carlo method and setting up the grid
   * @param rpt_fem_reconstruction_parameters Parameters for the particle position reconstruction
   * @param rpt_detector_parameters Parameters specifying information about the detectors for photon count calculation with the Monte Carlo method
   */
  RPTL2Projection(Parameters::RPTParameters &rpt_parameters,
                  Parameters::RPTFEMReconstructionParameters
                    &rpt_fem_reconstruction_parameters,
                  Parameters::DetectorParameters &rpt_detector_parameters,
                  Parameters::Mesh               &mesh_parameters)
    : mpi_communicator(MPI_COMM_WORLD)
    , fe(1)
    , dof_handler(triangulation)
    , mapping(FE_SimplexP<dim>(1))
    , rpt_parameters(rpt_parameters)
    , fem_reconstruction_parameters(rpt_fem_reconstruction_parameters)
    , detector_parameters(rpt_detector_parameters)
    , mesh_parameters(mesh_parameters)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)

  {}
  /**
   * @brief Using the L2 projection, builds a dictionary of photon counts
   * with respect to a detector for given positions inside a reactor.
   */
  void
  L2_project();

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
   * @brief Set up the particle_handler used for the data interpolation
   */
  void
  setup_particles();

  /**
   * @brief Assemble matrix and right hand side that form the linear system
   * that need to be solved for a given detector using the Beam Monte-Carlo
   * model.
   *
   * @param detector_no detector_no defines the detector for which the
   * linear system is being assembled.
   */
  void
  assemble_system_monte_carlo(unsigned detector_no);


  /**
   * @brief Assemble matrix and right hand side that form the linear system
   * that need to be solved for a given detector using raw experimental data.
   *
   * @param detector_no detector_no defines the detector for which the
   * linear system is being assembled.
   */
  void
  assemble_system_data(unsigned detector_no);

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
   * ".vtk" file format.
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
   * @brief Saves the built dictionary which will be used later for the 3D
   * FEM reconstruction.
   */
  void
  checkpoint();

  MPI_Comm                                   mpi_communicator;
  Triangulation<dim>                         triangulation;
  FE_SimplexP<dim>                           fe;
  DoFHandler<dim>                            dof_handler;
  MappingFE<dim>                             mapping;
  Parameters::RPTParameters                  rpt_parameters;
  Parameters::RPTFEMReconstructionParameters fem_reconstruction_parameters;
  Parameters::DetectorParameters             detector_parameters;
  Parameters::Mesh                           mesh_parameters;

  // In the data reconstruction mode, a particle handler is used to store the
  // information
  Particles::ParticleHandler<dim> particle_handler;


  std::vector<Detector<dim>> detectors;
  unsigned int               n_detector;

  AffineConstraints<double>   constraints;
  LA::MPI::SparseMatrix       system_matrix;
  LA::MPI::Vector             system_rhs;
  LA::MPI::Vector             solution;
  std::vector<Vector<double>> nodal_counts;

  ConditionalOStream pcout;
  TimerOutput        computing_timer;
};



template <int dim>
class RPTFEMReconstruction
{
public:
  /**
   * @brief Constructor for the RPTFEMReconstruction
   *
   * @param rpt_parameters General parameters for photon count calculation with the Monte Carlo method and setting up the grid
   * @param rpt_fem_reconstruction_parameters Parameters for the particle position reconstruction
   * @param rpt_detector_parameters Parameters specifying information about the detectors for photon count calculation with the Monte Carlo method
   */
  RPTFEMReconstruction(Parameters::RPTParameters &rpt_parameters,
                       Parameters::RPTFEMReconstructionParameters
                         &rpt_fem_reconstruction_parameters,
                       Parameters::DetectorParameters &rpt_detector_parameters,
                       Parameters::Mesh               &mesh_parameters)
    : fe(1)
    , dof_handler(triangulation)
    , mapping(FE_SimplexP<dim>(1))
    , parameters(rpt_parameters)
    , fem_reconstruction_parameters(rpt_fem_reconstruction_parameters)
    , detector_parameters(rpt_detector_parameters)
    , mesh_parameters(mesh_parameters)
    , computing_timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
  {}

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
   * @return cost It is calculated with the cost function mentioned in the
   * parameter file.
   */
  double
  calculate_cost(
    const TriaActiveIterator<DoFCellAccessor<dim, dim, false>> &cell,
    Vector<double>      &reference_location,
    const double        &last_constraint,
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
                                     const double   &last_constraint);

  /**
   * @brief Searches in the reference space the position of the particle.
   * If the found position respects the set tolerance, and the calculated
   * cost is minimal it calculates and stores the position of the particle
   * in system coordinates.
   *
   * @param count_from_all_detectors Contains the photon counts of all
   * detectors for all vertices of a cell.
   * @param experimental_count Experimental counts from all detectors for a
   * given particle position
   * @param cell Cell where the particle's position is being searched
   * @param tol_reference_location Tolerance on the validity of the found
   * location in the reference space
   * @param max_cost Cost that has to be minimized in order to find a
   * particle's real position
   * @param reference_location Position of the particle in reference
   * coordinates
   * @param position_found Bool,'true' if the particle's position was found,
   * and 'false' if the particle's position couldn't be found
   * @param real_location Position of the particle in system coordinates
   */
  void
  search_position_in_reference_space(
    std::vector<std::vector<double>> &count_from_all_detectors,
    std::vector<double>              &experimental_count,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const double   &tol_reference_location,
    double         &max_cost,
    Vector<double> &reference_location,
    bool           &position_found,
    Point<dim>     &real_location);

  /**
   * @brief Finds and saves the position of the particle in system coordinates.
   *
   * @param experimental_counts experimental_count contains the experimental
   * counts of every detector for a given position.
   *
   * @param tol_reference_location tolerance on the validity of the found
   * location in the reference space
   *
   * @return 'true' if the particle's position was found and 'false' if the
   * particle's position couldn't be found
   */
  bool
  find_position_global_search(std::vector<double> &experimental_counts,
                              const double         tol_reference_location);

  /**
   * @brief Finds and saves the position of the particle in system coordinates
   * by doing a local search around the previously found position's cell's
   * neighbors
   *
   * @param experimental_counts experimental_count contains the experimental
   * counts of every detector for a given position.
   *
   * @param tol_reference_location tolerance on the validity of the found
   * location in the reference space
   *
   * @return 'true' if the particle's position was found and 'false' if the
   * particle's position couldn't be found
   */
  bool
  find_position_local_search(std::vector<double> &experimental_counts,
                             const double         tol_reference_location);

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
   * @brief Loads the previously built dictionary for the 3D particle position
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
  // Work in progress to migrate the mesh functionalities
  // to use the regular mesh data container
  Parameters::Mesh mesh_parameters;



  std::vector<Detector<dim>>                     detectors;
  unsigned int                                   n_detector;
  std::vector<Vector<double>>                    nodal_counts;
  std::vector<Point<dim>>                        found_positions;
  typename DoFHandler<dim>::active_cell_iterator previous_position_cell;

  TimerOutput computing_timer;
};

/**
 * @brief Sets-up, assembles and solves the linear system to find the location
 * of a particle in reference coordinates.
 *
 * @param vertex_count contains the photon counts of all detectors all
 * vertices of a cell.
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
  std::vector<double>              &experimental_count,
  Parameters::RPTFEMReconstructionParameters::FEMCostFunction
    &cost_function_type);

#endif
