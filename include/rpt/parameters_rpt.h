/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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
 *  This file defines parameters needed for the RPT simulation in the Parameters
 * namespace.
 */

#ifndef lethe_rpt_parameters_h
#define lethe_rpt_parameters_h

#include <deal.II/base/config.h>

#include <core/parameters.h>

#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief RPTParameters - Defines the common parameters for the RPT
   * simulation as fixed parameters, number of Monte Carlo iterations and
   * filename of the particle positions.
   */

  struct RPTParameters
  {
    // Particle positions filename
    std::string particle_positions_file;

    // Show results in terminal during computation
    Parameters::Verbosity verbosity;

    // Enable to export counts result in a .csv file
    bool export_counts;

    // Exporting counts csv filename
    std::string export_counts_file;

    // Number of Monte Carlo iterations for alpha and theta
    unsigned int n_monte_carlo_iteration;

    // Seed of the random number generator
    int seed;

    // All parameters that are fixed by the user
    double reactor_radius; // [m]
    double reactor_height; // [m]
    double peak_to_total_ratio;
    double sampling_time; // [s]

    // Fixed parameters
    double gamma_rays_emitted; // Number of gamma-rays emitted by each
    // disintegration
    double attenuation_coefficient_detector; // Total linear attenuation
    // coefficient of the detector

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief InitialRPTParameters - Allows parameters tuning. If turned on,
   * other parameters values are the initial guess for the optimization,
   * otherwise those act as fixed parameters.
   */

  struct RPTTuningParameters
  {
    // Enable tuning parameters
    bool tuning;

    // Type of cost function
    enum class CostFunctionType
    {
      larachi,
      l1,
      l2
    };

    CostFunctionType cost_function_type;

    // Filename of experimental data
    std::string experimental_file;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief DetectorParameters - Defines information related to detectors. All
   * detectors must have the same dimensions (r and l).
   */

  struct DetectorParameters
  {
    double      radius;
    double      length;
    std::string detector_positions_file;

    // Tuned parameters
    std::vector<double>
      dead_time; // Dead time of the detector per accepted pulse
    std::vector<double> activity; // Activity of the tracer
    // Total linear attenuation coefficient of the medium
    std::vector<double> attenuation_coefficient_reactor;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct RPTReconstructionParameters
  {
    int reactor_refinement;
    int coarse_mesh_level; // Level of the coarse mesh where all the counts of
                           // the first vertices
    std::string reconstruction_counts_file;
    std::string reconstruction_positions_file;
    bool analyse_positions; // Allow to analyse results with known positions
    std::string known_positions_file;
    bool
      verbose_clock; // Allow to show total wallclock time elapsed since start

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief FEMReconstructionParameters - Defines the parameters for the
   * rpt_l2_projection_3d and rpt_fem_reconstruction_3d applications
   */
  struct RPTFEMReconstructionParameters
  {
    std::string  mesh_file;       // mesh filename
    unsigned int z_subdivisions;  // number of subdivisions of the initial grid
                                  // in z direction
    unsigned int mesh_refinement; // number of refinement the grid undergoes
    bool l2_project_and_reconstruct; // run the rpt_l2_projection_3d application
                                     // before reconstruction
    std::string experimental_counts_file; // file including experimental counts
                                          // from all detectors
    std::string export_positions_file;    // file including all found positions
    std::string dof_handler_file; // file with the saved DOFHandler object
    std::vector<std::string>
      nodal_counts_file; // vector containing the filenames of the files with
                         // the nodal counts from the built dictionary for each
                         // detector (1 file per detector)
    double extrapolation_tolerance; // tolerance when extrapolating from a cell
                                    // in the reference space to find the
                                    // particle's position


    std::string input_positions_file;
    std::string input_counts_file;

    // type of cost function applied when evaluating the particle's real
    // position
    enum class FEMCostFunction
    {
      absolute,
      relative
    };
    FEMCostFunction fem_cost_function;

    // type of mesh that is used to model the cylindrical vessel's geometry
    enum class FEMMeshType
    {
      dealii, // the grid used is a subdivided cylinder
      gmsh,
      dealiigen
    };
    FEMMeshType mesh_type;

    // type of searching algorithm
    enum class FEMSearchType
    {
      local,
      global
    };
    FEMSearchType search_type;

    // type of model used to build the L2 projection
    enum class FEMModel
    {
      monte_carlo,
      data
    };
    FEMModel model_type;

    unsigned int
         search_proximity_level; // level of proximity of the search scope
    bool verbose_clock_fem_reconstruction; // allow to show total wallclock time
                                           // elapsed since start

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

} // namespace Parameters



#endif // lethe_rpt_parameters_h
