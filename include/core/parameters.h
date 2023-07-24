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

#include <core/dimensionality.h>
#include <core/ib_particle.h>
#include <core/multiphysics.h>
#include <core/utilities.h>

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

  /** @brief Class to account for different fluid indicator:
   *  - fluid0: fluid 0 only,
   *  - fluid1: fluid 1 only,
   *  - both: both fluids
   * This is used for multiphase simulations, in PostProcessing and in
   * VOF (see parameter_multiphysics.h)
   */
  enum class FluidIndicator
  {
    fluid0,
    fluid1,
    both
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

    // Thermal conductivity liquid - Units in W/(m*K)
    double thermal_conductivity_l;

    // Thermal conductivity solid - Units in W/(m*K)
    double thermal_conductivity_s;

    // Thermal expansion coefficient liquid - Units in 1/K
    double thermal_expansion_l;

    // Thermal expansion coefficient solid - Units in 1/K
    double thermal_expansion_s;

    // kinematic viscosity liquid - Units in m^2/(s)
    double viscosity_l;

    // kinematic viscosity solid - Units in in m^2/(s)
    double viscosity_s;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality dimensions);
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
    // power-law does not allow for minimal viscosity
    double shear_rate_min;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality dimensions);
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
    parse_parameters(ParameterHandler &prm, const Dimensionality dimensions);
  };

  /**
   * @brief Non Newtonian - Defines the parameters for
   * non newtonian flows according to the chosen
   * rheological model.
   */
  struct NonNewtonian
  {
    CarreauParameters  carreau_parameters;
    PowerLawParameters powerlaw_parameters;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality dimensions);
  };


  /**
   * @brief Isothermal ideal gas model to solve for isothermal weakly compressible fluid
   * flows.
   */
  struct IsothermalIdealGasDensityParameters
  {
    // Reference state density of the gas in Pa
    double density_ref;
    // Specific gas constant in J/kg/K
    double R;
    // Absolute temperature of the ideal gas in K
    double T;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality dimensions);
  };

  /**
   * @brief SurfaceTensionParameters - Defines parameters for surface tension models
   */
  struct SurfaceTensionParameters
  {
    // Surface tension coefficient (sigma) in N/m
    double surface_tension_coefficient;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Material - Class that defines the physical property of a material.
   * Generally a material will be a fluid, but for conjugated heat transfer,
   * this may also be a solid.
   */
  class Material
  {
  public:
    Material()
    {}

    void
    declare_parameters(ParameterHandler &prm,
                       std::string       material_prefix,
                       unsigned int      id);
    void
    parse_parameters(ParameterHandler &   prm,
                     std::string          material_prefix,
                     const unsigned int   id,
                     const Dimensionality dimensions);

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
    // tracer diffusivity in L^2/s
    double tracer_diffusivity;

    // Phase change parameters
    PhaseChange phase_change_parameters;

    // Non Newtonian model parameters
    enum class RheologicalModel
    {
      powerlaw,
      carreau,
      newtonian,
      phase_change
    } rheological_model;
    NonNewtonian non_newtonian_parameters;

    enum class DensityModel
    {
      constant,
      isothermal_ideal_gas
    } density_model;
    IsothermalIdealGasDensityParameters isothermal_ideal_gas_density_parameters;

    enum class SpecificHeatModel
    {
      constant,
      phase_change
    } specific_heat_model;

    enum class ThermalConductivityModel
    {
      constant,
      linear,
      phase_change
    } thermal_conductivity_model;

    enum class ThermalExpansionModel
    {
      constant,
      phase_change
    } thermal_expansion_model;

    // Linear thermal conductivity parameters: k = k_A0 + k_A1 * T
    double k_A0;
    double k_A1;
  };

  /**
   * @brief MaterialInteractions - Class that defines physical properties due to interactions between two different materials (either fluid-fluid or fluid-solid).
   */
  class MaterialInteractions
  {
  public:
    MaterialInteractions()
    {}

    enum class MaterialInteractionsType
    {
      fluid_fluid,
      fluid_solid
    } material_interaction_type;

    // Surface tension models
    enum class SurfaceTensionModel
    {
      constant
    } surface_tension_model;
    SurfaceTensionParameters surface_tension_parameters;

    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_fluid_interaction_with_material_interaction_id;
    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_solid_interaction_with_material_interaction_id;

    void
    declare_parameters(ParameterHandler &prm, unsigned int id);

    void
    parse_parameters(ParameterHandler &prm, const unsigned int id);
  };

  /**
   * @brief PhysicalProperties - Define the possible physical properties.
   * All continuum equations share the same physical properties object but only
   * take the subset of properties they require
   * Defined as a class with public attributes in order to use a non-static
   * declare_parameters methods (useful for multiple fluid simulations).
   */
  class PhysicalProperties
  {
  public:
    PhysicalProperties()
    {}

    // Fluid objects for multiphase simulations
    std::vector<Material>     fluids;
    unsigned int              number_of_fluids;
    static const unsigned int max_fluids = 2;

    // Solid objects for conjugated simulations
    std::vector<Material>     solids;
    unsigned int              number_of_solids;
    static const unsigned int max_solids = 1;

    // Fluid-fluid or fluid-solid interactions
    std::vector<MaterialInteractions> material_interactions;
    unsigned int                      number_of_material_interactions;
    static const unsigned int         max_material_interactions = 3;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_fluid_interactions_with_material_interaction_ids;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_solid_interactions_with_material_interaction_ids;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &   prm,
                     const Dimensionality dimensions = Dimensionality());
  };

  /**
   * @brief Stabilization - Defines parameters for an advanced control over the stabilization strategy used by the solvers.
   */
  struct Stabilization
  {
    bool use_default_stabilization;

    // pressure scaling factor used to facilitate the linear solving when
    // pressure and velocity have very different scales
    double pressure_scaling_factor;

    enum class NavierStokesStabilization
    {
      pspg_supg,
      gls,
      grad_div
    } stabilization;

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
   * @brief Laser_FreeSurfaceRadiation - Defines the subparameters for
   * the radiation sink at the free surface (air/metal interface).
   */
  struct Laser_FreeSurfaceRadiation
  {
    bool enable_radiation;

    // Parameters for the radiation term at the melt pool free surface
    double Stefan_Boltzmann_constant;
    double emissivity;
    double Tinf;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Laser parameters - Defines the parameters for the
   * laser heat source.
   */
  template <int dim>
  class Laser
  {
  public:
    // A boolean parameter that enables the calculations of laser heat source
    bool activate_laser;

    // Laser concentration factor indicates the definition of the beam radius.
    // In almost all the articles, it is assumed equal to 2.0
    double concentration_factor;

    // Laser power in W
    double laser_power;

    // Absorptivity is defined as the fraction of the amount of incident
    // radiation that is absorbed by the surface, and it is measured using
    // diffuse reﬂectance spectroscopy (DRS). Generally, a constant value in
    // the range of 0.3-0.8 (for welding processes with titanium) are used
    // in the literature. However, recent studies show that it varies with
    // powder particle size distribution, and the angle of incidence that
    // changes due to the dynamic melt pool surface [Zhang, Z., Huang, Y.,
    // Kasinathan, A.R., Shahabad, S.I., Ali, U., Mahmoodkhani, Y. and
    // Toyserkani, E., 2019. 3-Dimensional heat transfer modeling for laser
    // powder-bed fusion additive manufacturing with volumetric heat sources
    // based on varied thermal conductivity and absorptivity. Optics & Laser
    // Technology, 109, pp.297-312.]
    double laser_absorptivity;

    // Penetration depth of the laser
    double penetration_depth;

    // Laser beam radius on the melt pool surface
    double beam_radius;

    // Beam orientation shows the orientation of the laser beam. For instance,
    // if a laser beam is emitted perpendicular on a plane in x-y coordinates,
    // the orientation of the laser beam will be in the z direction. Note that
    // this parameter cannot be equal to z in two-dimensional simulations.
    // Plus and minus shows the direction of the laser beam
    enum class BeamOrientation
    {
      x_plus,
      y_plus,
      z_plus,
      x_minus,
      y_minus,
      z_minus,
    } beam_orientation;

    // beam_orientation_coordinate parameter stores the integer (x = 0, y = 1,
    // z =2) value of the beam_orientation parameter
    unsigned int beam_orientation_coordinate;

    // beam_direction shows the direction of laser beam (either in positive (1)
    // or negative (0) direction
    bool beam_direction;

    // Based on the laser beam orientation, the integer values of a
    // perpendicular plane to the laser beam orientation are stored in the
    // following parameters (x = 0, y = 1, z = 2)
    unsigned int perpendicular_plane_coordinate_one;
    unsigned int perpendicular_plane_coordinate_two;

    // Laser scan path indicates the path of the laser focal point during a
    // simulation
    std::shared_ptr<Functions::ParsedFunction<dim>> laser_scan_path;

    // Start and end time of the laser operation
    double start_time;
    double end_time;

    // Parameters for the radiation term at the melt pool free surface
    Laser_FreeSurfaceRadiation radiation;

    void
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

    // Enable flow rate post-processing
    bool calculate_flow_rate;

    // Set initial time to start calculations for velocities
    double initial_time;

    // Frequency of the calculation of the post-processed quantity
    unsigned int calculation_frequency;

    // Frequency of the output
    unsigned int output_frequency;

    // Prefix for kinetic energy output
    std::string kinetic_energy_output_name;

    // Prefix for pressure drop output
    std::string pressure_drop_output_name;

    // Prefix for flow rate output
    std::string flow_rate_output_name;

    // Prefix for the enstrophy output
    std::string enstrophy_output_name;

    // Prefix for the apparent viscosity output
    std::string apparent_viscosity_output_name;

    // Enable tracer statistics
    bool calculate_tracer_statistics;

    // Prefix for the tracer output
    std::string tracer_output_name;

    // Enable temperature statistics
    bool calculate_temperature_statistics;

    // Prefix for the temperature output
    std::string temperature_output_name;

    // Enable heat flux calculation
    bool calculate_heat_flux;

    // Prefix for the total heat flux output
    std::string heat_flux_output_name;

    // Fluid domain, used when post-processing a multiphase simulation
    Parameters::FluidIndicator postprocessed_fluid;

    // Enable barycenter calculation for fluid 1 in VOF simulations
    bool calculate_vof_barycenter;

    // Prefix for the barycenter output
    std::string barycenter_output_name;

    // Enable smoothing postprocessed vectors and scalars
    bool smoothed_output_fields;

    // Enable phase statistics
    bool calculate_phase_statistics;

    // Prefix for the phase output
    std::string phase_output_name;

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

    // Interpolation order Cahn-Hilliard
    unsigned int phase_ch_order;
    unsigned int potential_ch_order;

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

    // Abort solver if non-linear solution has not reached tolerance
    bool abort_at_convergence_failure;


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
      periodic_hills,
      cylinder
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

    // A boolean parameter which enables adding the neighbor boundary cells of
    // boundary cells in DEM simulations. This parameter should only be enabled
    // for simulations with concave geometries (for instance particles inside a
    // drum). In simulations with convex geometries, it must not be enabled.
    // This is also reported to users in a warning in
    // find_boundary_cells_information.
    bool expand_particle_wall_contact_search;

    // Grid displacement at initiation
    Tensor<1, 3> translation;

    // Grid rotation at initiation
    Tensor<1, 3> rotation_axis;
    double       rotation_angle;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Mesh Adaptation Parameters that can differ for each refinement variable
   */
  struct MultipleAdaptationParameters
  {
    // Coarsening fraction
    double coarsening_fraction;

    // Refinement fraction
    double refinement_fraction;
  };

  /**
   * @brief MeshAdaption - Parameters that control dynamic mesh adaptation.
   * Dynamic mesh adaptation in Lethe is very flexible and can be both local
   * and global.
   */
  struct MeshAdaptation
  {
    // Initial adaptive refinement
    unsigned int initial_refinement;

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
      phase,
      temperature,
      phase_ch,
      chemical_potential_ch
    } variable;

    // Map containing the refinement variables
    std::map<Variable, MultipleAdaptationParameters> variables;
    // declaration for parsing variables
    Variable                     vars;
    MultipleAdaptationParameters var_adaptation_param;

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

    // Enable the control of the mesh refinement to target a specific number of
    // elements equal to the maximum number of elements.
    bool mesh_controller_is_enabled;

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

    enum class TestType
    {
      particles,
      mobility_status,
      subdomain
    } test_type;

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
    unsigned int                 levels_not_precalculated;
    double                       inside_radius;
    double                       outside_radius;
    bool                         time_extrapolation_of_refinement_zone;
    std::vector<IBParticle<dim>> particles;
    bool                         calculate_force_ib;
    bool                         assemble_navier_stokes_inside;
    std::string                  ib_force_output_file;
    Tensor<1, dim>               gravity;
    double                       particle_nonlinear_tol;
    double                       wall_youngs_modulus;
    double                       wall_poisson_ratio;
    double                       wall_rolling_friction_coefficient;
    double                       wall_friction_coefficient;
    double                       wall_restitution_coefficient;
    unsigned int                 coupling_frequency;
    bool                         enable_lubrication_force;
    double                       lubrication_range_max;
    double                       lubrication_range_min;
    double                       contact_search_radius_factor;
    int                          contact_search_frequency;
    bool                         load_particles_from_file;
    std::string                  particles_file;
    bool                         enable_extra_sharp_interface_vtu_output_field;

    std::shared_ptr<Functions::ParsedFunction<dim>> f_gravity;

    double      particle_nonlinear_tolerance;
    double      length_ratio;
    double      alpha;
    bool        print_dem;
    std::string ib_particles_pvd_file;
  };

  /**
   * @brief FlowControl - Set average velocity on a boundary (CFD) or the domain
   * (CFD-DEM).
   */
  struct DynamicFlowControl
  {
    // Enable flow control
    bool enable_flow_control;

    // Average velocity target (L/t)
    double average_velocity_0;

    // Boundary id at flow inlet
    unsigned int boundary_flow_id;

    // Flow direction (x=0, y=1 ,z=2)
    unsigned int flow_direction;

    // Initial beta
    double beta_0;

    // Relaxation coefficient for beta force controller
    // beta_n+1 = beta_n + alpha * (...)
    double alpha;

    // If beta at n+1 step is in this threshold over beta at n step, beta n+1
    // is kept as beta n. This avoid a new term of force in the matrix
    double beta_threshold;

    // Type of verbosity for the flow control
    Verbosity verbosity;

    // Apply scaled beta force to particles
    bool enable_beta_particle;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

} // namespace Parameters
#endif
