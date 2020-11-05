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
 * Author: Simon Gauvin, Bruno Blais, Polytechnique Montreal, 2019-
 */

/*
 * This file defines the parameter namespace. This namespace
 * contains the classical structures which are used to structure
 * the various simulations that can be carried out using Lethe.
 * The parameters structures are constructed in logical building
 * blocks so that a solver can adequately choose which blocks
 * are required.
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <core/ib_particle.h>

#ifndef lethe_parameters_h
#  define lethe_parameters_h

using namespace dealii;

namespace Parameters
{
  enum class Verbosity
  {
    quiet,
    verbose
  };

  /**
   * @brief SimulationControl - Defines the parameter that control the flow of the simulation
   * as well as the frequency of the output of the solutions.
   */

  struct SimulationControl
  {
    // Method used for time progression of eulerian solvers (steady, unsteady)
    enum class TimeSteppingMethod
    {
      steady,
      steady_bdf,
      bdf1,
      bdf2,
      bdf3,
      sdirk2,
      sdirk2_1,
      sdirk2_2,
      sdirk3,
      sdirk3_1,
      sdirk3_2,
      sdirk3_3
    } method;

    // Method used for time progression (steady, unsteady)
    enum class LagrangianTimeSteppingMethod
    {
      explicit_euler,
      velocity_verlet
    } lagrangian_method;

    // Initial time step
    double dt;

    // End time
    double timeEnd;

    // Adaptative time stepping
    bool adapt;

    // Max CFL
    double maxCFL;

    // Aimed tolerance at which simulation is stopped
    double stop_tolerance;

    // Max CFL
    double adaptative_time_step_scaling;

    // BDF startup time scaling
    double startup_timestep_scaling;

    // Number of mesh adaptation (steady simulations)
    unsigned int number_mesh_adaptation;

    // Folder for simulation output
    std::string output_folder;

    // Prefix for simulation output
    std::string output_name;

    enum class OutputControl
    {
      iteration,
      time
    } output_control;

    // Frequency of the output
    unsigned int output_frequency;

    // Frequency of the output
    double output_time;

    // Enable output of the boundaries
    bool output_boundaries;

    // Frequency of the log output to the terminal
    unsigned int log_frequency;

    // Display precision of the log output to the terminal
    unsigned int log_precision;

    // Subdivisions of the results in the output
    unsigned int subdivision;

    // Subdivisions of the results in the output
    unsigned int group_files;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief PhysicalProperties - Define the possible physical properties.
   * All continuum equations share the same physical properties object but only
   * take the subset of properties they require
   */
  struct PhysicalProperties
  {
    // Kinematic viscosity (mu/rho) in units of L^2/s
    double viscosity;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief Timer - Defines the parameters that control the timing of the simulation.
   * Lethe supports advanced timing features which supports the monitoring of
   * specific sub-sections of the software to evaluate the relative
   * workload.
   */
  struct Timer
  {
    // Time measurement in the simulation. None, at each iteration, only at the
    // end
    enum class Type
    {
      none,
      iteration,
      end
    };

    Type type;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Forces - Defines the parameters for the
   * force calculation on boundaries of the domain.
   */
  struct Forces
  {
    // Type of verbosity for the force calculation
    Verbosity verbosity;

    // Enable force post-processing
    bool calculate_force;

    // Enable torque post-processing
    bool calculate_torque;

    // Frequency of the output
    unsigned int calculation_frequency;

    // Frequency of the output
    unsigned int output_frequency;

    // Output precision
    unsigned int output_precision;

    // Prefix for simulation output
    std::string force_output_name;

    // Prefix for the torque output
    std::string torque_output_name;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Postprocessing - Defines the parameters
   * for the post-processing. In Lethe, post-processing
   * implies the calculation of quantities derived from the principal
   * variables. For example, the integral kinetic energy or the integral
   * enstrophy.
   */
  struct PostProcessing
  {
    Verbosity verbosity;

    // Enable total kinetic energy post-processing
    bool calculate_kinetic_energy;

    // Enable total enstrophy post-processing
    bool calculate_enstrophy;

    // Enable velocity post-processing
    bool calculate_average_velocities;

    // Enable nondimensionalize average velocities
    bool nondimensionalization;

    // Set initial time to start calculations for velocities
    double initial_time;

    // Id of the boundary where the flow inlet
    unsigned int id_flow_control;

    // Flow direction (x=0, y=1 ,z=2)
    unsigned int flow_direction;

    // Component value of all the time-averaged velocities (often called
    // simply average_velocities) is averaged at the component_location
    double component_average;

    // Component value where the time-averaged velocities (often called
    // simply average_velocities) is averaged
    double component_location;

    // Frequency of the calculation of the post-processed quantity
    unsigned int calculation_frequency;

    // Frequency of the output
    unsigned int output_frequency;

    // Prefix for kinectic energy output
    std::string kinetic_energy_output_name;

    // Prefix for the enstrophy output
    std::string enstrophy_output_name;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief FEM - The finite element section
   * controls the properties of the finite element method. This section
   * constrols the order of polynomial integration and the number of quadrature
   * points within the cells.
   */
  struct FEM
  {
    // Interpolation order velocity
    unsigned int velocity_order;

    // Interpolation order pressure
    unsigned int pressure_order;

    // Number of quadrature points per dimension
    // The final number of quadrature point will be
    // number_quadrature_points^dim
    unsigned int number_quadrature_points;

    // Apply high order mapping everywhere
    bool qmapping_all;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief NonLinearSolver - Parameter that controls the solution of the
   * non-linear problems.
   */
  struct NonLinearSolver
  {
    // Type of linear solver
    enum class SolverType
    {
      newton,
      skip_newton,
      adaptative_newton
    };

    Verbosity verbosity;

    // Type of non-linear solver
    SolverType solver;

    // Tolerance
    double tolerance;

    // Relative Tolerance
    double step_tolerance;

    // Maximal number of iterations for the Newton solver
    unsigned int max_iterations;

    // Residual precision
    unsigned int display_precision;

    // Iterations to skip in the non-linear solver
    unsigned int skip_iterations;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief LinearSolver - Parameters that controls the solution of the
   * linear system of equations that arise from the finite element problem.
   */
  struct LinearSolver
  {
    // Type of linear solver
    enum class SolverType
    {
      gmres,
      bicgstab,
      amg,
      direct
    };
    SolverType solver;

    Verbosity verbosity;

    // Relative residual of the iterative solver
    double relative_residual;

    // Minimum residual of the iterative solver
    double minimum_residual;

    // Maximum number of iterations
    int max_iterations;

    // Maximum number of krylov vectors
    int max_krylov_vectors;

    // ILU or ILUT fill
    double ilu_precond_fill;

    // ILU or ILUT absolute tolerance
    double ilu_precond_atol;

    // ILU or ILUT relative tolerance
    double ilu_precond_rtol;

    // ILU or ILUT fill
    double amg_precond_ilu_fill;

    // ILU or ILUT absolute tolerance
    double amg_precond_ilu_atol;

    // ILU or ILUT relative tolerance
    double amg_precond_ilu_rtol;

    // AMG aggregation threshold
    double amg_aggregation_threshold;

    // AMG number of cycles
    unsigned int amg_n_cycles;

    // AMG W_cycle
    bool amg_w_cycles;

    // AMG Smoother sweeps
    unsigned int amg_smoother_sweeps;

    // AMG Smoother overalp
    unsigned int amg_smoother_overlap;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Mesh - Parameters that control mesh reading and mesh generation.
   */
  struct Mesh
  {
    // GMSH or dealii
    enum class Type
    {
      gmsh,
      dealii,
      periodic_hills
    };
    Type type;

    // File name of the mesh
    std::string file_name;

    // Name of the grid in GridTools
    std::string grid_type;

    // Arguments of the GridTools
    std::string grid_arguments;

    // Initial refinement level of primitive mesh
    unsigned int initial_refinement;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief MeshAdaption - Parameters that control dynamic mesh adaptation.
   * Dynamic mesh adaptation in Lethe is very flexible and can be both local
   * and global.
   */
  struct MeshAdaptation
  {
    // Type of mesh adaptation
    enum class Type
    {
      none,
      uniform,
      kelly
    } type;

    enum class Variable
    {
      velocity,
      pressure
    } variable;

    // Decision factor for Kelly refinement (number or fraction)
    enum class FractionType
    {
      number,
      fraction
    } fractionType;

    // Maximum number of elements
    unsigned int maximum_number_elements;

    // Maximum refinement level
    unsigned int maximum_refinement_level;

    // Minimum refinement level
    unsigned int minimum_refinement_level;

    // Refinement after frequency iter
    unsigned int frequency;

    // Refinement fraction
    double refinement_fraction;

    // Coarsening fraction
    double coarsening_fraction;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief Testing - Some solvers have a specific testing
   * mode that can be enabled to output more variables to the
   * terminal. This is enabled using the Testing parameters.
   */
  struct Testing
  {
    bool enabled;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief Restart - Controls writing and reading
   * simulation checkpoints.
   */

  struct Restart
  {
    std::string  filename;
    bool         restart;
    bool         checkpoint;
    unsigned int frequency;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief VelocitySource - Adds velocity-dependent
   * source term to the Navier-Stokes equations with
   * the appropriate jacobian matrix. Currently only
   * a change to a rotating frame is supported, but additional
   * terms such as a Darcy force or similar could be easily
   * added.
   */

  struct VelocitySource
  {
    enum class VelocitySourceType
    {
      none,
      srf
    };
    VelocitySourceType type;
    double             omega_x;
    double             omega_y;
    double             omega_z;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  template <int dim>
  class IBParticles
  {
  public:
    unsigned int                 nb;
    unsigned int                 order;
    unsigned int                 nb_force_eval;
    unsigned int                 initial_refinement;
    double                       inside_radius;
    double                       outside_radius;
    std::vector<IBParticle<dim>> particles;
    bool                         calculate_force_ib;
    std::string                  ib_force_output_file;


    static void
    declare_parameters(ParameterHandler &prm);
    static void
    declare_default_entry(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief FlowControl - Set volumetric flow rate on a boundary id
   * toward the normal direction of this wall.
   */
  struct DynamicFlowControl
  {
    // Enable flow control
    bool enable_flow_control;

    // Volumetric flow rate (L^3/t)
    double flow_rate;

    // Id of the boundary where the flow inlet
    unsigned int id_flow_control;

    // Flow direction (x=0, y=1 ,z=2)
    unsigned int flow_direction;

    // Initial beta
    double beta_0;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


} // namespace Parameters
#endif
