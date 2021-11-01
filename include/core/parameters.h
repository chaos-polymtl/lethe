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



#ifndef lethe_parameters_h
#define lethe_parameters_h

#include <core/ib_particle.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>



using namespace dealii;

namespace Parameters
{
  enum class Verbosity
  {
    quiet,
    verbose,
    extra_verbose
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
      sdirk22,
      sdirk22_1,
      sdirk22_2,
      sdirk33,
      sdirk33_1,
      sdirk33_2,
      sdirk33_3
    } method;

    // Method used for time progression (steady, unsteady)
    enum class LagrangianTimeSteppingMethod
    {
      explicit_euler,
      velocity_verlet,
      gear3
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

    enum class BDFStartupMethods
    {
      initial_solution,
      sdirk_step,
      multiple_step_bdf,
    } bdf_startup_method;


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
   * @brief Fluid - Class for fluid definition
   */
  class Fluid
  {
  public:
    Fluid()
    {}

    void
    declare_parameters(ParameterHandler &prm, unsigned int id);
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    // Kinematic viscosity (nu = mu/rho) in units of L^2/s
    double viscosity;
    // volumetric mass density (rho) in units of kg/m^3
    double density;
    // specific heat capacity (cp) in J/K/kg
    double specific_heat;
    // thermal conductivity (k) in W/m/K
    double thermal_conductivity;
    // thermal expansion coefficient (alpha) in 1/K
    double thermal_expansion;
    // tracer diffusivity) in L^2/s
    double tracer_diffusivity;
  };

  /**
   * @brief Power-law rheological model to solve for non Newtonian
   * flows.
   */
  struct PowerLawParameters
  {
    // Fluid consistency index
    double K;
    // Flow behavior index"
    double n;
    // Minimal shear rate magnitude for which we calculate viscosity, since
    // power-law does not allow for minimal visocsity
    double shear_rate_min;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Carreau rheological model to solve for non Newtonian
   * flows.
   */
  struct CarreauParameters
  {
    // Viscosity of the flow when the shear rate tends to 0
    double viscosity_0;
    // Hypothetical viscosity of the flow when the shear rate is very high
    double viscosity_inf;
    // Relaxation time
    double lambda;
    // Carreau parameter
    double a;
    // Power parameter
    double n;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Non Newtonian - Defines the parameters for
   * non newtonian flows according to the chosen
   * rheological model.
   */

  struct NonNewtonian
  {
    // Non Newtonian model
    enum class Model
    {
      powerlaw,
      carreau
    } model;

    CarreauParameters  carreau_parameters;
    PowerLawParameters powerlaw_parameters;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief Phase change model for melting/freezing liquids
   * The model assumes that the phase change occurs between
   * a solidus and liquidus temperature. This defines a solidification
   * interval which is used to smooth the non-linearity of the melting problem
   * or to fit the real thermodynamics of the melting process.
   *
   */
  struct PhaseChange
  {
    // Solidus temperature - Units in K
    double T_solidus;

    // Liquidus temperature - Units in K
    double T_liquidus;

    // Latent enthalpy for the phase change - Units in J/kg
    double latent_enthalpy;

    // Specific heat liquid - Units in J/(kg*K)
    double cp_l;

    // Specific heat solid - Units in J/(kg*K)
    double cp_s;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief InterfaceSharpening - Defines the parameters for
   * interface sharpening in the VOF solver.
   */
  struct InterfaceSharpening
  {
    // Interface sharpening parameters. The sharpening method and parameters are
    // explained in the dam break VOF example:
    // https://github.com/lethe-cfd/lethe/wiki/Dam-break-VOF
    // sharpening_threshold is the phase fraction threshold for sharpening. It
    // should be chosen in the range of (0,1), but generally it is equal to 0.5
    // interface_sharpness is a parameter which defines the sharpness of the
    // interface. It should be chosen in the range of (1,2] sharpening_frequency
    // (integer) is the frequency at which the interface sharpneing is called.
    // Users may set this variable to 1 to call interface sharpening at every
    // step, but it could be chosen in the range of [1-20]

    double sharpening_threshold;
    double interface_sharpness;
    int    sharpening_frequency;
    // Type of verbosity for the interface sharpening calculation
    Verbosity verbosity;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief PhysicalProperties - Define the possible physical properties.
   * All continuum equations share the same physical properties object but only
   * take the subset of properties they require
   * Defined as a class with public attributes in order to use a non-static
   * declare_paremeters methods (useful for multiple fluid simulations).
   */
  class PhysicalProperties
  {
  public:
    PhysicalProperties()
    {}

    // Non Newtonian parameters
    bool         non_newtonian_flow;
    NonNewtonian non_newtonian_parameters;

    // Phase change parameters
    bool        enable_phase_change;
    PhaseChange phase_change_parameters;

    // Fluid objects for multiphasic simulations
    std::vector<Fluid>        fluids;
    unsigned int              number_of_fluids;
    static const unsigned int max_fluids = 2;

    void
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

    bool write_time_in_error_table;

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

    // Enable calculating apparent viscosity
    bool calculate_apparent_viscosity;

    // Enable velocity post-processing
    bool calculate_average_velocities;

    // Enable pressure drop post-processing
    bool calculate_pressure_drop;

    // The outlet boundary ID for pressure drop calculation
    unsigned int inlet_boundary_id;

    // The outlet boundary ID for pressure drop calculation
    unsigned int outlet_boundary_id;

    // Set initial time to start calculations for velocities
    double initial_time;

    // Frequency of the calculation of the post-processed quantity
    unsigned int calculation_frequency;

    // Frequency of the output
    unsigned int output_frequency;

    // Prefix for kinectic energy output
    std::string kinetic_energy_output_name;

    // Prefix for pressure drop output
    std::string pressure_drop_output_name;

    // Prefix for the enstrophy output
    std::string enstrophy_output_name;

    // Prefix for the apparent viscosity output
    std::string apparent_viscosity_output_name;

    // Enable tracer statistics
    bool calculate_tracer_statistics;

    // Prefix for the tracer output
    std::string tracer_output_name;

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

    // Interpolation order void fraction
    unsigned int void_fraction_order;

    // Interpolation order temperature
    unsigned int temperature_order;

    // Interpolation order tracer
    unsigned int tracer_order;

    // Interpolation order vof model
    unsigned int VOF_order;

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
    // Type of non-linear solver
    enum class SolverType
    {
      newton,
      inexact_newton,
      kinsol_newton
    };

    // Kinsol solver strategy
    enum class KinsolStrategy
    {
      normal_newton,
      line_search,
      picard
    };

    Verbosity verbosity;

    // Type of non-linear solver
    SolverType solver;

    // Kinsol solver strategy
    KinsolStrategy kinsol_strategy;

    // Tolerance
    double tolerance;

    // Maximal number of iterations for the Newton solver
    unsigned int max_iterations;

    // Residual precision
    unsigned int display_precision;


    // Force RHS recalculation at the beginning of every non-linear steps
    // This is required if there is a fixed point component to the non-linear
    // solver that is changed at the beginning of every newton iteration.
    // This is notably the case of the sharp edge method.
    // The default value of this parameter is false.
    bool force_rhs_calculation;

    // Matrix reconstruction tolerance
    // This parameter controls the reconstruction of the system matrix
    // If the residual after a newton step is lower than previous_residual *
    // matrix_tolerance then that iteration is considered sufficient and the
    // matrix is not reassembled at the next iteration.
    double matrix_tolerance;

    // Relative Tolerance
    double step_tolerance;

    // Carry jacobian matrix over to the new non-linear problem
    bool reuse_matrix;


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

    // Block linear solver to throw error.
    bool force_linear_solver_continuation;

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

    // Enabling fixing initial refinement from a target size
    bool refine_until_target_size;

    // Allowing the use of a simplex mesh
    bool simplex;

    // Target size when automatically refining initial mesh
    double target_size;

    // Enables checking the input grid for diamond-shaped cells
    bool check_for_diamond_cells;

    // Grid displacement at initiation
    bool   translate;
    double delta_x;
    double delta_y;
    double delta_z;

    // Grid rotation at initiation
    bool   rotate;
    int    axis;
    double angle;

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
      pressure,
      phase
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

  struct MeshBoxRefinement
  {
    // GMSH or dealii
    std::shared_ptr<Mesh> box_mesh;
    // Initial refinement level of primitive mesh contained in the box
    unsigned int initial_refinement;

    void
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
    void
    declare_parameters(ParameterHandler &prm);
    void
    declare_default_entry(ParameterHandler &prm, unsigned int index);
    void
    parse_parameters(ParameterHandler &prm);

    unsigned int                 nb;
    unsigned int                 order;
    unsigned int                 initial_refinement;
    double                       inside_radius;
    double                       outside_radius;
    std::vector<IBParticle<dim>> particles;
    bool                         calculate_force_ib;
    bool                         assemble_navier_stokes_inside;
    std::string                  ib_force_output_file;
    double                       density;
    Tensor<1, dim>               gravity;
    double                       particle_nonlinear_tol;
    double                       wall_youngs_modulus;
    double                       wall_poisson_ratio;
    double                       wall_rolling_friction_coefficient;
    double                       wall_friction_coefficient;
    double                       wall_restitution_coefficient;
    unsigned int                 coupling_frequency;

    std::shared_ptr<Functions::ParsedFunction<dim>> f_gravity;


    double      particle_nonlinear_tolerance;
    double      length_ratio;
    double      alpha;
    bool        integrate_motion;
    std::string ib_particles_pvd_file;
  };

  /**
   * @brief FlowControl - Set volumetric flow rate on a boundary id
   * toward the normal direction of this wall.
   */
  struct DynamicFlowControl
  {
    // Type of verbosity for the flow control
    Verbosity verbosity;

    // Enable flow control
    bool enable_flow_control;

    // Volumetric flow rate (L^3/t)
    double flow_rate_0;

    // Boundary id at flow inlet
    unsigned int boundary_flow_id;

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
