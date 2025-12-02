// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

using namespace dealii;

namespace Parameters
{
  struct SizeOfSubsections
  {
    int boundary_conditions;
    int manifolds;
  };


  /**
   * @brief Extract the maximum number of all variable size sections within the parameter file
   *
   * @param file_name Name of the parameter file from which the size are parsed
   */
  SizeOfSubsections
  get_size_of_subsections(const std::string &file_name);

  enum class Verbosity : std::uint8_t
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
  enum class FluidIndicator : std::uint8_t
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
    enum class TimeSteppingMethod : std::uint8_t
    {
      steady,
      steady_bdf,
      bdf1,
      bdf2,
      bdf3,
      sdirk22,
      sdirk33,
      sdirk43
    } method;

    // Method used for time progression (steady, unsteady)
    enum class LagrangianTimeSteppingMethod : std::uint8_t
    {
      explicit_euler,
      velocity_verlet,
      gear3
    } lagrangian_method;

    // Initial time step
    double dt;

    // End time
    double time_end;

    // Boolean to keep the time step for the last iteration regardless of the
    // end time specify. Both for fixed time step and adaptive time step.
    bool time_step_independent_of_end_time;

    // Adaptative time stepping
    bool adapt;

    // Max CFL
    double maxCFL;

    // Max time step
    double max_dt;

    // Aimed tolerance at which simulation is stopped
    double stop_tolerance;

    // Rate of increase of the time step value
    double adaptative_time_step_scaling;

    // BDF startup time scaling
    double startup_timestep_scaling;

    // True if the time-step should be overriden upon restart
    bool override_time_step_on_restart;

    // Number of mesh adaptation (steady simulations)
    unsigned int number_mesh_adaptation;

    // Folder for simulation output
    std::string output_folder;

    // Prefix for simulation output
    std::string output_name;

    enum class OutputControl : std::uint8_t
    {
      iteration,
      time
    } output_control;

    enum class BDFStartupMethods : std::uint8_t
    {
      initial_solution,
      multiple_step_bdf,
    } bdf_startup_method;

    // Frequency of the output (for iteration output control)
    unsigned int output_iteration_frequency;

    // Time frequency of the output (for time output control)
    double output_time_frequency;

    // Output at specific times (for time output control)
    std::vector<double> output_times_vector;

    // Time window for file output (for both iteration and time output control)
    std::vector<double> output_time_interval;

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

    // Specific heat of liquid - Units in J/(kg*K)
    double cp_l;

    // Specific heat of solid - Units in J/(kg*K)
    double cp_s;

    // Thermal conductivity of liquid - Units in W/(m*K)
    double thermal_conductivity_l;

    // Thermal conductivity of solid - Units in W/(m*K)
    double thermal_conductivity_s;

    // Thermal expansion coefficient of liquid - Units in 1/K
    double thermal_expansion_l;

    // Thermal expansion coefficient of solid - Units in 1/K
    double thermal_expansion_s;

    // kinematic viscosity of liquid - Units in m^2/(s)
    double kinematic_viscosity_l;

    // kinematic viscosity of solid - Units in m^2/(s)
    double kinematic_viscosity_s;

    // Darcy penalty of liquid - Units in 1/(s)
    double penalty_l;

    // Darcy penalty of solid - Units in 1/(s)
    double penalty_s;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Power-law rheological model to solve for non-Newtonian
   * flows.
   */
  struct PowerLawParameters
  {
    // Fluid consistency index
    double K;
    // Flow behavior index
    double n;
    // Minimal shear rate magnitude for which we calculate kinematic viscosity,
    // since power-law does not allow for minimal kinematic viscosity
    double shear_rate_min;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Carreau model to solve for non-Newtonian
   * flows.
   */
  struct CarreauParameters
  {
    // Kinematic viscosity of the flow when the shear rate tends to 0
    double kinematic_viscosity_0;
    // Hypothetical kinematic viscosity of the flow when the shear rate is very
    // high
    double kinematic_viscosity_inf;
    // Relaxation time
    double lambda;
    // Carreau parameter
    double a;
    // Power parameter
    double n;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief non-Newtonian - Defines the parameters for non-Newtonian flows
   * according to the chosen rheological model.
   */
  struct NonNewtonian
  {
    CarreauParameters  carreau_parameters;
    PowerLawParameters powerlaw_parameters;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Tanh-based physical properties model to handle properties in/out of immersed solids
   */
  struct ImmersedSolidTanhParameters
  {
    // Properties that apply to the tracer physics with immersed solids
    double tracer_diffusivity_inside;
    double tracer_diffusivity_outside;
    double tracer_reaction_constant_inside;
    double tracer_reaction_constant_outside;
    double thickness;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Gaussian-based physical properties model to handle properties at the
   * interface of immersed solids and in the bulk of each phase
   */
  struct ImmersedSolidGaussianParameters
  {
    // Properties that apply to the tracer physics with immersed solids
    double tracer_diffusivity_interface;
    double tracer_diffusivity_bulk;
    double tracer_reaction_constant_interface;
    double tracer_reaction_constant_bulk;
    double thickness;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Isothermal ideal gas model to solve for isothermal weakly
   * compressible fluid flows.
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
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief SurfaceTensionParameters - Defines parameters for surface tension
   * models
   */
  struct SurfaceTensionParameters
  {
    // Surface tension coefficient (sigma or sigma_0) in N/m
    double surface_tension_coefficient;
    // Temperature of the reference state corresponding to the surface tension
    // coefficient (T_0) in K
    double T_0;
    // Surface tension gradient with respect to the temperature (dsigma/dT) in
    // N/(m*K)
    double surface_tension_gradient;

    // Solidus temperature - Units in K
    double T_solidus;

    // Liquidus temperature - Units in K
    double T_liquidus;

    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    static void
    declare_parameters(ParameterHandler &prm);
    void
    /**
     * @brief Parse the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     *
     * @param[in] dimensions The Dimensionality object controlling the
     * fundamental dimensions (length, time, mass, temperature) of the problem.
     */
    parse_parameters(const ParameterHandler           &prm,
                     const Parameters::Dimensionality &dimensions);
  };

  /**
   * @brief MobilityCahnHilliardParameters - Defines parameters for the mobility
   * models used in the Cahn-Hilliard equations.
   */
  struct MobilityCahnHilliardParameters
  {
    // Mobility constant (M) in m^2/s
    double mobility_cahn_hilliard_constant;

    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     *
     * @param[in] dimensions The Dimensionality object controlling the
     * fundamental dimensions (length, time, mass, temperature) of the problem.
     */
    void
    parse_parameters(const ParameterHandler           &prm,
                     const Parameters::Dimensionality &dimensions);
  };


  /**
   * @brief Material - Class that defines the physical property of a material.
   * Generally a material will be a fluid, but for conjugated heat transfer,
   * this may also be a solid.
   */
  class Material
  {
  public:
    Material() = default;

    void
    declare_parameters(ParameterHandler  &prm,
                       const std::string &material_prefix,
                       unsigned int       id) const;
    void
    parse_parameters(ParameterHandler     &prm,
                     const std::string    &material_prefix,
                     const unsigned int    id,
                     const Dimensionality &dimensions);

    // Kinematic viscosity (nu = mu/rho) in units of L^2/s
    double kinematic_viscosity;
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

    // tracer reaction constant in 1/s^[order]
    double tracer_reaction_constant;
    // tracer reaction order
    double tracer_reaction_order;
    // When the reaction order <1, concentrations are used in the denominator of
    // a few terms. This threshold prevents the introduction of NaN entries as
    // long as it is set above 0.
    double tracer_reaction_threshold;

    // Phase change parameters
    PhaseChange phase_change_parameters;

    // non-Newtonian model parameters
    enum class RheologicalModel : std::int8_t
    {
      powerlaw,
      carreau,
      newtonian,
      phase_change
    } rheological_model;
    NonNewtonian non_newtonian_parameters;

    enum class DensityModel : std::int8_t
    {
      constant,
      isothermal_ideal_gas
    } density_model;
    IsothermalIdealGasDensityParameters isothermal_ideal_gas_density_parameters;

    enum class SpecificHeatModel : std::int8_t
    {
      constant,
      phase_change
    } specific_heat_model;

    enum class ThermalConductivityModel : std::int8_t
    {
      constant,
      linear,
      phase_change
    } thermal_conductivity_model;

    enum class ThermalExpansionModel : std::int8_t
    {
      constant,
      phase_change
    } thermal_expansion_model;

    enum class TracerDiffusivityModel : std::int8_t
    {
      constant,
      immersed_boundary_tanh,
      immersed_boundary_gaussian
    } tracer_diffusivity_model;

    enum class TracerReactionPrefactorModel : std::int8_t
    {
      none,
      constant,
      immersed_boundary_tanh,
      immersed_boundary_gaussian
    } tracer_reaction_prefactor_model;

    // Struct that contains the parameters to handle physical properties when
    // immersed solids tanh models are used
    ImmersedSolidTanhParameters immersed_solid_tanh_parameters;

    // Struct that contains the parameters to handle physical properties when
    // immersed solids Gaussian models are used
    ImmersedSolidGaussianParameters immersed_solid_gaussian_parameters;

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
    MaterialInteractions() = default;

    enum class MaterialInteractionsType : std::int8_t
    {
      fluid_fluid,
      fluid_solid
    } material_interaction_type;

    // Surface tension models
    enum class SurfaceTensionModel : std::int8_t
    {
      constant,
      linear,
      phase_change
    } surface_tension_model;
    SurfaceTensionParameters surface_tension_parameters;

    // Cahn-Hilliard mobility models
    enum class MobilityCahnHilliardModel : std::int8_t
    {
      constant,
      quartic
    } mobility_cahn_hilliard_model;
    MobilityCahnHilliardParameters mobility_cahn_hilliard_parameters;

    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_fluid_interaction_with_material_interaction_id;
    std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_solid_interaction_with_material_interaction_id;

    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     *
     * @param[in] id The material id.
     */
    void
    declare_parameters(ParameterHandler &prm, unsigned int id) const;

    /**
     * @brief Parse the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     *
     * @param[in] id The material id.
     *
     * @param[in] dimensions The Dimensionality object controlling the
     * fundamental dimensions (length, time, mass, temperature) of the problem.
     */
    void
    parse_parameters(ParameterHandler                 &prm,
                     const unsigned int                id,
                     const Parameters::Dimensionality &dimensions);
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
    PhysicalProperties() = default;

    // Fluid objects for multiphase simulations
    std::vector<Material>         fluids;
    unsigned int                  number_of_fluids;
    static constexpr unsigned int max_fluids = 2;

    // Solid objects for conjugated simulations
    std::vector<Material>         solids;
    unsigned int                  number_of_solids;
    static constexpr unsigned int max_solids = 1;

    // Fluid-fluid or fluid-solid interactions
    std::vector<MaterialInteractions> material_interactions;
    unsigned int                      number_of_material_interactions;
    static constexpr unsigned int     max_material_interactions = 3;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_fluid_interactions_with_material_interaction_ids;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_solid_interactions_with_material_interaction_ids;

    /*
     * Reference Temperature for all physical properties of fluids and solids.
     * Currently, this is only used by the thermal expansion models.
     */
    double reference_temperature;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler     &prm,
                     const Dimensionality &dimensions = Dimensionality());
  };

  /**
   * @brief Set of parameters constraining a certain portion of a fluid domain
   * to null velocity and pressure fields to mimic a solid subdomain.
   *
   * @remark Pressure DOFs in "solid" cells that are next to "fluid" cells are
   * not constrained.
   *
   * @tparam dim Number of dimensions of the problem (2D or 3D).
   *
   * @note At the moment, only the temperature field is used to constrain the
   * "solid" domain.
   */
  template <int dim>
  struct ConstrainSolidDomain
  {
    /// Enable/disable (@p true/false) the solid domain constraining feature.
    bool enable;

    /// Total number of constraints (maximum of 1 per fluid)
    unsigned int number_of_constraints;

    /// Identifiers of fluids that are constrained
    std::vector<unsigned int> fluid_ids;

    /// Absolute tolerance applied on filtered phase fraction for
    /// constrained cell selection
    std::vector<double> filtered_phase_fraction_tolerance;

    /// Lower threshold values of the constraining field (temperature)
    std::vector<double> temperature_min_values;

    /// Upper threshold values of the constraining field (temperature)
    std::vector<double> temperature_max_values;

    /// Enable/disable (@p true/false) the definition of a plane for geometrical
    /// restrictions on the domain where the stasis constraint is applied.
    bool enable_domain_restriction_with_plane;

    /// Coordinates of a point on the restriction plane for the stasis
    /// constraint application domain
    Point<dim> restriction_plane_point;

    /// Outward-pointing normal vector from the restricted domain to define the
    /// restriction plane for the stasis constraint application domain
    Tensor<1, dim> restriction_plane_normal_vector;

    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm ParameterHandler object.
     *
     * @param[in] max_number_of_constraints Maximum number of zero velocity
     * constraints applied to the domain.
     */
    void
    declare_parameters(ParameterHandler  &prm,
                       const unsigned int max_number_of_constraints);

    /**
     * @brief Parse the parameters.
     *
     * @param[in] prm ParameterHandler object.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /**
     *
     * @brief Declare the default parameters for each constraint.
     *
     * @param[in,out] prm ParameterHandler object.
     *
     */
    static void
    declare_default_entries(ParameterHandler &prm);

    /**
     *
     * @brief Parse parameters for each constraint.
     *
     * @param[in] prm ParameterHandler object.
     *
     * @param[in] constraint_id Identifiers of the constraint (1 per fluid). The
     * numbering starts at 0.
     */
    void
    parse_constraint_parameters(const ParameterHandler &prm,
                                const unsigned int      constraint_id);
  };

  /**
   * @brief Stabilization - Defines parameters for an advanced control over the stabilization strategy used by the solvers.
   */
  struct Stabilization
  {
    // Defines if default stabilization parameters should be used
    bool use_default_stabilization;

    bool heat_transfer_dcdd_stabilization;

    /// Boolean indicating if the DCDD stabilization term for the VOF phase
    /// fraction should be assembled (@p true) or not (@p false).
    bool vof_dcdd_stabilization;

    // Pressure scaling factor used to facilitate the linear solving when
    // pressure and velocity have very different scales
    double pressure_scaling_factor;

    enum class NavierStokesStabilization : std::int8_t
    {
      pspg_supg,
      gls,
      grad_div
    } stabilization;

    enum class ScalarLimiters : std::int8_t
    {
      moe,
      none
    } scalar_limiter;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Timer - Defines the parameters that control the timing of the simulation.
   * Lethe supports advanced timing features which supports the monitoring of
   * specific subsections of the software to evaluate the relative
   * workload.
   */
  struct Timer
  {
    // Time measurement in the simulation. None, at each iteration, only at the
    // end
    enum class Type : std::int8_t
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

    // Type of laser model used in simulations. With "exponential_decay", the
    // laser acts as a volumetric source, whereas, with
    // "heat_flux_vof_interface", the laser behaves as a surface flux at the
    // interface between fluids (VOF auxiliary physic must be enabled to use
    // this model).
    enum class LaserType : std::int8_t
    {
      exponential_decay,
      gaussian_heat_flux_vof_interface,
      uniform_heat_flux_vof_interface
    } laser_type;

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
    enum class BeamOrientation : std::int8_t
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

    // beam_direction shows the direction of laser beam (either in positive
    // (true) or negative (false) direction
    bool beam_direction;

    // Based on the laser beam orientation, the integer values of a
    // perpendicular plane to the laser beam orientation are stored in the
    // following parameters (x = 0, y = 1, z = 2)
    unsigned int perpendicular_plane_coordinate_one;

    // Beam axis
    Tensor<1, dim> beam_axis;

    unsigned int perpendicular_plane_coordinate_two;

    // rotation angle of the laser axis in rad
    double rotation_angle;

    // rotation axis
    Tensor<1, dim> rotation_axis;

    // rotation matrix
    Tensor<2, dim> rotation_matrix;

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

    /// Enable total kinetic energy post-processing
    bool calculate_kinetic_energy;

    /// Enable total enstrophy post-processing
    bool calculate_enstrophy;

    /// Enable pressure power post-processing
    bool calculate_pressure_power;

    /// Enable viscous dissipation post-processing
    bool calculate_viscous_dissipation;

    /// Enable calculating apparent viscosity
    bool calculate_apparent_viscosity;

    /// Enable velocity post-processing
    bool calculate_average_velocities;

    /// Enable average temperature and average heat flux post-processing
    bool calculate_average_temp_and_hf;

    /// Enable pressure drop post-processing
    bool calculate_pressure_drop;

    /// The inlet boundary ID for pressure drop calculation
    unsigned int inlet_boundary_id;

    /// The outlet boundary ID for pressure drop calculation
    unsigned int outlet_boundary_id;

    /// Enable flow rate post-processing
    bool calculate_flow_rate;

    /// Enable tracer flow rate post-processing
    bool calculate_tracer_flow_rate;

    /// Set initial time to start calculations for velocities
    double initial_time_for_average_velocities;

    /// Set initial time to start calculations for average temperature and
    /// average heat flux
    double initial_time_for_average_temp_and_hf;

    /// Frequency of the calculation of the post-processed quantity
    unsigned int calculation_frequency;

    /// Frequency of the output
    unsigned int output_frequency;

    /// Prefix for kinetic energy output
    std::string kinetic_energy_output_name;

    /// Prefix for pressure drop output
    std::string pressure_drop_output_name;

    /// Prefix for flow rate output
    std::string flow_rate_output_name;

    /// Prefix for tracer flow rate output
    std::string tracer_flow_rate_output_name;

    /// Prefix for the enstrophy output
    std::string enstrophy_output_name;

    /// Prefix for the pressure power output
    std::string pressure_power_output_name;

    /// Prefix for the viscous dissipation output
    std::string viscous_dissipation_output_name;

    /// Prefix for the apparent viscosity output
    std::string apparent_viscosity_output_name;

    /// Enable tracer statistics
    bool calculate_tracer_statistics;

    /// Prefix for the tracer output
    std::string tracer_output_name;

    /// Enable temperature statistics
    bool calculate_temperature_statistics;

    /// Prefix for the temperature output
    std::string temperature_output_name;

    /// Enable calculation of liquid fraction in phase change problems
    bool calculate_liquid_fraction;

    /// Prefix for the temperature output
    std::string liquid_fraction_output_name;

    /// Enable heat flux calculation
    bool calculate_heat_flux;

    /// Prefix for the total heat flux output
    std::string heat_flux_output_name;

    /// Fluid domain, used when post-processing a multiphase simulation
    Parameters::FluidIndicator postprocessed_fluid;

    /// Enable barycenter calculation for fluid 1 in VOF and Cahn-Hilliard
    /// simulations
    bool calculate_barycenter;

    /// Prefix for the VOF and Cahn-Hilliard barycenter output
    std::string barycenter_output_name;

    /// Enable phase statistics
    bool calculate_phase_statistics;

    /// Prefix for the phase output
    std::string phase_output_name;

    /// Enable mass conservation calculation for both fluids in VOF simulations
    bool calculate_mass_conservation;

    /// Prefix for the VOF mass conservation output
    std::string mass_conservation_output_name;

    /// Enable energies calculation on the domain in Cahn-Hilliard simulations
    bool calculate_phase_energy;

    /// Prefix for the energy output in Cahn-Hilliard simulations
    std::string phase_energy_output_name;

    /// Enable calculation of total fluid volume and total particles volume in
    /// cfd-dem simulation
    bool calculate_phase_volumes;

    /// prefix for the total volume output in cfd-dem simulation
    std::string phase_volumes_output_name;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief FEM - The finite element section
   * controls the properties of the finite element method. This section
   * controls the order of polynomial integration and the number of quadrature
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

    // Switch tracer to DG formulation instead of CG
    bool tracer_uses_dg;

    // Interpolation order vof model
    unsigned int VOF_order;

    // Switch vof to DG formulation instead of CG
    bool VOF_uses_dg;

    // Interpolation order Cahn-Hilliard
    unsigned int phase_cahn_hilliard_order;
    unsigned int potential_cahn_hilliard_order;

    // Option for bubble enrichment functions
    bool enable_bubble_function_velocity;
    bool enable_bubble_function_pressure;

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
    enum class SolverType : std::int8_t
    {
      newton,
      inexact_newton,
      kinsol_newton
    };

    // Kinsol solver strategy
    enum class KinsolStrategy : std::int8_t
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

    // Reuse preconditioner for the next non-linear iterations
    bool reuse_preconditioner;

    // Abort solver if non-linear solution has not reached tolerance
    bool abort_at_convergence_failure;

    static void
    declare_parameters(ParameterHandler &prm, const std::string &physics_name);
    void
    parse_parameters(ParameterHandler &prm, const std::string &physics_name);
  };

  /**
   * @brief LinearSolver - Parameters that control the solution of the
   * linear system of equations that arise from the finite element problem for
   * each of the physics available in Lethe
   */
  struct LinearSolver
  {
    // Type of linear solver
    enum class SolverType : std::int8_t
    {
      gmres,
      bicgstab,
      direct
    };

    SolverType solver;

    /// Verbosity of linear solver
    Verbosity verbosity;

    /// Flag to rescale linear and non-linear residuals by the sqrt of the
    /// triangulation volume
    bool rescale_residual_by_volume;

    /// Relative residuals of the iterative solver
    double relative_residual;

    /// Minimum residual of the iterative solver
    double minimum_residual;

    /// Maximum number of iterations
    int max_iterations;

    /// Maximum number of krylov vectors
    int max_krylov_vectors;

    /// Enable hessians in jacobian
    bool enable_hessians_jacobian;

    /// Enable hessians in residual
    bool enable_hessians_residual;

    /// Type of preconditioner
    enum class PreconditionerType : std::int8_t
    {
      ilu,
      amg,
      lsmg,
      gcmg
    };
    PreconditionerType preconditioner;

    /// ILU or ILUT fill
    unsigned int ilu_precond_fill;

    /// ILU or ILUT absolute tolerance
    double ilu_precond_atol;

    /// ILU or ILUT relative tolerance
    double ilu_precond_rtol;

    /// AMG parameters either as linear solver preconditioner or as
    /// preconditioner of a coarse-grid solver for LSMG or GCMG

    /// ILU or ILUT fill for smoother
    unsigned int amg_precond_ilu_fill;

    /// ILU or ILUT absolute tolerance for smoother
    double amg_precond_ilu_atol;

    /// ILU or ILUT relative tolerance for smoother
    double amg_precond_ilu_rtol;

    /// AMG aggregation threshold
    double amg_aggregation_threshold;

    /// AMG number of cycles
    unsigned int amg_n_cycles;

    /// AMG W_cycle
    bool amg_w_cycles;

    /// AMG Smoother sweeps
    unsigned int amg_smoother_sweeps;

    /// AMG Smoother overlap
    unsigned int amg_smoother_overlap;

    /// Block linear solver to throw error.
    bool force_linear_solver_continuation;

    /// MG min level
    int mg_min_level;

    /// MG minimum number of cells per level
    int mg_level_min_cells;

    /// MG intermediate level
    int mg_int_level;

    /// MG enable hessians in jacobian
    bool mg_enable_hessians_jacobian;

    /// Type of multigrid
    enum class MultigridCoarseningSequenceType : std::int8_t
    {
      h,
      p,
      hp,
      ph
    };
    MultigridCoarseningSequenceType mg_coarsening_type;

    /// Type of p coarsening sequence
    MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
      mg_p_coarsening_type;

    /// MG smoother number of iterations
    int mg_smoother_iterations;

    /// MG smoother relaxation parameter
    double mg_smoother_relaxation;

    /// Type of preconditioner for the MG smoother
    enum class MultigridSmootherPreconditionerType : std::int8_t
    {
      InverseDiagonal,
      AdditiveSchwarzMethod
    };
    MultigridSmootherPreconditionerType mg_smoother_preconditioner_type;

    /// MG eigenvalue estimation for smoother relaxation parameter
    bool mg_smoother_eig_estimation;

    /// MG smoothing range to set range between eigenvalues
    double eig_estimation_smoothing_range;

    /// MG number of cg iterations to find eigenvalue
    int eig_estimation_cg_n_iterations;

    /// MG print max, min, eigenvalues
    Verbosity eig_estimation_verbose;

    /// Type of coarse grid solver
    enum class CoarseGridSolverType : std::int8_t
    {
      gmres,
      amg,
      ilu,
      direct
    };
    CoarseGridSolverType mg_coarse_grid_solver;

    /// MG use FE_Q_iso_Q1 elements for coarse grid
    bool mg_use_fe_q_iso_q1;

    /// MG coarse-grid solver maximum number of iterations
    int mg_gmres_max_iterations;

    /// MG coarse-grid solver tolerance
    double mg_gmres_tolerance;

    /// MG coarse-grid solver reduce
    double mg_gmres_reduce;

    /// MG coarse-grid solver maximum number of krylov vectors
    int mg_gmres_max_krylov_vectors;

    /// MG coarse-grid solver preconditioner
    PreconditionerType mg_gmres_preconditioner;

    /// MG use default parameters for AMG
    bool mg_amg_use_default_parameters;

    /// MG information about levels
    Verbosity mg_verbosity;

    static void
    declare_parameters(ParameterHandler &prm, const std::string &physics_name);
    void
    parse_parameters(ParameterHandler &prm, const std::string &physics_name);
  };

  /**
   * @brief Mesh - Parameters that control mesh reading and mesh generation.
   */
  struct Mesh
  {
    // GMSH or dealii
    enum class Type : std::int8_t
    {
      gmsh,
      dealii,
      periodic_hills,
      cylinder
    };
    Type type;

    /// File name of the mesh
    std::string file_name;

    /// Name of the grid in GridTools
    std::string grid_type;

    /// Arguments of the GridTools
    std::string grid_arguments;

    /// Initial refinement level of primitive mesh
    unsigned int initial_refinement;

    /// Initial refinement level of primitive mesh near user-defined boundary
    /// conditions
    unsigned int initial_refinement_at_boundaries;

    /// List of boundary ids to refine
    std::vector<int> boundaries_to_refine;

    /// Enable fixing initial refinement from a target size
    bool refine_until_target_size;

    /// Allow the use of a simplex mesh
    bool simplex;

    /// Target size when automatically refining initial mesh
    double target_size;

    /// Enable checking the input grid for diamond-shaped cells
    bool check_for_diamond_cells;

    /* A boolean parameter which enables adding the neighbor boundary cells of
    * boundary cells in DEM simulations. This parameter should only be enabled
    * for simulations with concave geometries (for instance particles inside a
    * drum). In simulations with convex geometries, it must not be enabled.
    * This is also reported to users in a warning in
     find_boundary_cells_information.*/
    bool expand_particle_wall_contact_search;

    // Grid displacement at initiation
    Tensor<1, 3> translation;

    // Grid rotation at initiation
    Tensor<1, 3> rotation_axis;
    double       rotation_angle;

    /// Rescale the grid by the scale factor
    double scale;

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
    enum class Type : std::int8_t
    {
      none,
      uniform,
      kelly
    } type;

    /// Fields on which the mesh adaptation can be based
    Variable variable;

    // Map containing the refinement variables
    std::map<Variable, MultipleAdaptationParameters> variables;
    // declaration for parsing variables
    Variable                     vars;
    MultipleAdaptationParameters var_adaptation_param;

    // Decision factor for Kelly refinement (number or fraction)
    enum class FractionType : std::int8_t
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
    bool         refinement_at_frequency;

    // Enable the control of the mesh refinement to target a specific number of
    // elements equal to the maximum number of elements.
    bool mesh_controller_is_enabled;

    // Specifies if mesh adaptation can be used on certain boundaries
    bool is_boundary_refinement_fixed;

    // List of boundary ids to fix to their initial refinement state
    std::vector<int> boundaries_to_fix;

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

    enum class TestType : std::int8_t
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
    enum class DarcySourceType : std::int8_t
    {
      none,        // No Darcy source term
      phase_change // Phase change darcy source term which applies a
                   // penalization depending on the phase change model
    };

    enum class RotatingFrameType : std::int8_t
    {
      none,
      srf
    };

    RotatingFrameType rotating_frame_type;
    double            omega_x;
    double            omega_y;
    double            omega_z;

    /*
     * Type of darcy velocity source term applied to the Navier-Stokes equations
     */
    DarcySourceType darcy_type;


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
    // Vector of particles
    std::vector<IBParticle<dim>> particles;

    // Number of declared IB particles
    unsigned int nb_particles;
    // Boolean to determine whether the Navier-Stokes equations are
    // solved inside the particles.
    bool assemble_navier_stokes_inside;

    // Polynomial order of the IB stencil
    unsigned int order;
    // The length ratio used for the stencil calculation of the IB condition.
    double length_ratio;
    // Boolean controlling whether extrapolation is used to impose the
    // immersed boundary condition. If false, the IB condition is directly
    // imposed using nearest neighbors. the immersed boundary condition or not.
    // If it is set to false, all cut cells are fully imposed on the IB.
    bool enable_extrapolation;

    // Boolean for the calculation of the force at the IB
    bool calculate_force_ib;
    // Boolean for extra vtu field output
    bool enable_extra_sharp_interface_vtu_output_field;
    // The name of the output file for the forces on the IB.
    std::string ib_force_output_file;
    // Particles pvd file name
    std::string ib_particles_pvd_file;
    // Boolean for printing DEM information
    bool print_dem;

    // Number of initial refinements around each particle
    unsigned int initial_refinement;
    // Inner radius factor of the refinement zone
    double inside_radius;
    // Outer radius factor of the refinement zone
    double outside_radius;
    // Boolean for time-dependent mesh refinement according to the particle's
    // current position the refinement zone.
    bool time_extrapolation_of_refinement_zone;

    // Number of DEM time steps per CFD time step.
    unsigned int coupling_frequency;
    // Relaxation parameter for the CFD-DEM coupling.
    double alpha;
    // Frequency at which the contact search is performed at the CFD time
    // scale (once every X CFD time steps)
    int contact_search_frequency;
    // Particles' radius multiplier used to calculate the effective radius of
    // contact search
    double contact_search_radius_factor;
    // Boolean for lubrication force
    bool enable_lubrication_force;
    // The maximal range for which the lubrication force is evaluated. This
    // variable multiplies the smallest cell diameter to obtain the actual
    // range.
    double lubrication_range_max;
    // The minimal range for which the lubrication force is evaluated. This
    // variable multiplies the smallest cell diameter to obtain the actual
    // range.
    double lubrication_range_min;
    // Tolerance for the particle dynamics nonlinear solver.
    double particle_nonlinear_tolerance;
    // Boolean for the explicit calculation of the contact impulse in the
    // CFD-DEM coupling impulsion in the CFD-DEM coupling.
    bool explicit_contact_impulsion_calculation;
    // Boolean for explicit evaluation of the particle's position in the CFD-DEM
    // coupling
    bool explicit_position_integration_calculation;
    // Boolean for approximation of the contact radius. If true, the effective
    // radius replaces the actual particle's local curvature radius in the
    // calculation.
    bool approximate_radius_for_contact;
    // Function defining the gravitational acceleration vector used by the
    // CFD-DEM calculation.
    std::shared_ptr<Functions::ParsedFunction<dim>> f_gravity;
    // Young's modulus of the wall
    double wall_youngs_modulus;
    // Poisson ratio of the wall
    double wall_poisson_ratio;
    // Rolling friction coefficient of the wall
    double wall_rolling_friction_coefficient;
    // Sliding friction coefficient of the wall
    double wall_friction_coefficient;
    // Coefficient of restitution of the wall
    double wall_restitution_coefficient;

    // Boolean for loading particles from an independent file. If true, the
    // definition of the particles in the Particle subsection in the .prm file
    // is ignored. of the particle subsection.
    bool load_particles_from_file;
    // Name of the independent file containing particles' information at
    // insertion. Only used if load_particles_from_file is true.
    std::string particles_file;
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
    // is kept as beta n. This avoids a new term of force in the matrix
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

  /**
   * @brief Evaporation - Defines the subparameters for
   * the evaporation cooling and recoil pressure at the free surface
   * (air/metal interface).
   */
  struct Evaporation
  {
    enum class EvaporativeMassFluxModelType : std::int8_t
    {
      constant,
      temperature_dependent
    } evaporative_mass_flux_model_type;

    bool enable_evaporation_cooling;
    bool enable_recoil_pressure;

    // Parameters for the evaporation terms at the melt pool free surface
    double evaporation_mass_flux;
    double evaporation_coefficient;
    double recoil_pressure_coefficient;
    double molar_mass;
    double boiling_temperature;
    double latent_heat_evaporation;
    double ambient_pressure;
    double ambient_gas_density;
    double liquid_density;
    double universal_gas_constant;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Mortar - Defines parameters used to construct mortar elements.
   */
  template <int dim>
  struct Mortar
  {
    /// Indicates whether mortar elements are enabled
    bool enable;
    /// Mesh parameters for the rotor part
    std::shared_ptr<Mesh> rotor_mesh;
    /// Boundary ID # of the rotor at the rotor-stator interface
    unsigned int rotor_boundary_id;
    /// Boundary ID # of the stator at the rotor-stator interface
    unsigned int stator_boundary_id;
    /// Center of rotation of the rotor domain
    Point<dim> center_of_rotation;
    /// Rotation axis of the rotor domain
    Tensor<1, dim> rotation_axis;
    /// Rotation angle of the rotor domain in radians
    std::shared_ptr<Functions::ParsedFunction<dim>> rotor_rotation_angle;
    /// Angular velocity of the rotor domain
    std::shared_ptr<Functions::ParsedFunction<dim>> rotor_angular_velocity;
    /// Penalty factor for mortar elements
    double sip_factor;
    /// Oversampling factor for quadrature points
    unsigned int oversampling_factor;
    /// Type of verbosity for mortar
    Verbosity verbosity;

    void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
  };


} // namespace Parameters
#endif
