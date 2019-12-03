/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Simon Gauvin, Polytechnique Montreal, 2019
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#ifndef LETHE_GLS_PARAMETERS_H
#  define LETHE_GLS_PARAMETERS_H

using namespace dealii;

namespace Parameters
{
  enum Verbosity
  {
    quiet,
    verbose
  };

  struct SimulationControl
  {
    // Method used for time progression (steady, unsteady)
    enum TimeSteppingMethod
    {
      steady,
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
    };
    TimeSteppingMethod method;

    // Initial time step
    double dt;

    // End time
    double timeEnd;

    // Adaptative time stepping
    bool adapt;

    // Max CFL
    double maxCFL;

    // BDF startup time scaling
    double startup_timestep_scaling;

    // Number of mesh adaptation (steady simulations)
    unsigned int nbMeshAdapt;

    // Folder for simulation output
    std::string output_folder;

    // Prefix for simulation output
    std::string output_name;

    // Frequency of the output
    unsigned int outputFrequency;

    // Subdivisions of the results in the output
    unsigned int subdivision;

    // Subdivisions of the results in the output
    unsigned int group_files;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct PhysicalProperties
  {
    // Kinematic viscosity (mu/rho)
    double viscosity;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Timer
  {
    // Time measurement in the simulation. None, at each iteration, only at the
    // end
    enum Type
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

  struct Forces
  {
    // Type of verbosity for the iterative solver

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

    // Display precision
    unsigned int display_precision;

    // Prefix for simulation output
    std::string force_output_name;

    // Prefix for the torque output
    std::string torque_output_name;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  struct PostProcessing
  {
    Verbosity verbosity;

    // Enable total kinetic energy post-processing
    bool calculate_kinetic_energy;

    // Enable total enstrophy post-processing
    bool calculate_enstrophy;

    // Frequency of the output
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

  struct FEM
  {
    // Interpolation order velocity
    unsigned int velocityOrder;

    // Interpolation order pressure
    unsigned int pressureOrder;

    // Number of quadrature points
    unsigned int quadraturePoints;

    // Apply high order mapping everywhere
    bool qmapping_all;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct NonLinearSolver
  {
    // Type of linear solver
    enum SolverType
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

  struct LinearSolver
  {
    // Type of linear solver
    enum SolverType
    {
      gmres,
      bicgstab,
      amg
    };
    SolverType solver;

    Verbosity verbosity;

    // Residual precision
    unsigned int residual_precision;

    // Relative residual of the iterative solver
    double relative_residual;

    // Minimum residual of the iterative solver
    double minimum_residual;

    // Maximum number of iterations
    int max_iterations;

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

  struct Mesh
  {
    // GMSH or dealii primitive
    enum Type
    {
      gmsh,
      dealii,
      primitive
    };
    Type type;

    // Primitive types
    enum PrimitiveType
    {
      hyper_cube,
      hyper_shell,
      cylinder
    };
    PrimitiveType primitiveType;

    bool colorize;

    double arg1;
    double arg2;
    double arg3;
    double arg4;
    double arg5;
    double arg6;

    // File name of the mesh
    std::string file_name;

    // Name of the grid in GridTools
    std::string grid_type;

    // Arguments of the GridTools
    std::string grid_arguments;

    // Initial refinement level of primitive mesh
    unsigned int initialRefinement;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct MeshAdaptation
  {
    // Type of mesh adaptation
    enum Type
    {
      none,
      uniform,
      kelly
    };
    Type type;

    enum Variable
    {
      velocity,
      pressure
    };
    Variable variable;

    // Decision factor for KELLY refinement (number or fraction)
    enum FractionType
    {
      number,
      fraction
    };
    FractionType fractionType;

    // Maximum number of elements
    unsigned int maxNbElements;

    // Maximum refinement level
    unsigned int maxRefLevel;

    // Maximum refinement level
    unsigned int minRefLevel;

    // Refinement after frequency iter
    unsigned int frequency;

    // Refinement fractioni havent used ILUT much)
    double fractionRefinement;

    // Coarsening fraction
    double fractionCoarsening;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Testing
  {
    // Time measurement in the simulation. None, at each iteration, only at the
    // end
    bool enabled;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct Restart
  {
    // Time measurement in the simulation. None, at each iteration, only at the
    // end
    std::string  filename;
    bool         restart;
    bool         checkpoint;
    unsigned int frequency;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

} // namespace Parameters
#endif
