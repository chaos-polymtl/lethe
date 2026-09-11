// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
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

#include <limits>

using namespace dealii;

namespace Parameters
{
  struct SizeOfSubsections
  {
    int boundary_conditions = 0;
    int manifolds           = 0;
  };


  /**
   * @brief Extract the maximum number of all variable size sections within the parameter file
   *
   * @param[in] file_name Name of the parameter file from which the sizes are
   * parsed
   *
   * @param[in] require_subsection_size If true, the parameter file must
   * indicate the size of at least one variable size subsection, which is the
   * case for every solver for which a "boundary conditions" subsection is
   * mandatory. The DEM solvers, whose boundary conditions are declared in a
   * "DEM boundary conditions" subsection with a "number of boundary
   * conditions" entry, set this to false and consequently obtain a size of
   * zero for parameter files that do not declare any manifold.
   *
   * //TODO make the size of the DEM boundary conditions parsed instead.
   *
   * @return The size of the variable size subsections of the parameter file.
   */
  SizeOfSubsections
  get_size_of_subsections(const std::string &file_name,
                          const bool         require_subsection_size = true);

  enum class Verbosity : std::uint8_t
  {
    quiet,
    verbose,
    extra_verbose
  };

  /**
   * @brief Class to account for different fluid indicators.
   * This is used for multiphase simulations in PostProcessing, in
   * CLS (see parameter_multiphysics.h), and in VelocitySource (phase change
   * Carman-Kozeny permeability model).
   */
  enum class FluidIndicator : std::uint8_t
  {
    fluid0, ///< fluid 0 only
    fluid1, ///< fluid 1 only
    both    ///< both fluids
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
    } method = TimeSteppingMethod::steady;

    // Method used for time progression (steady, unsteady)
    enum class LagrangianTimeSteppingMethod : std::uint8_t
    {
      explicit_euler,
      velocity_verlet,
      gear3
    } lagrangian_method = LagrangianTimeSteppingMethod::explicit_euler;

    // Initial time step
    double dt = 1.;

    // End time
    double time_end = 1.;

    // End iteration (only used when end_control is EndControl::iteration)
    unsigned int iteration_end = 10;

    // Boolean to keep the time step for the last iteration regardless of the
    // end time specify. Both for fixed time step and adaptive time step.
    bool time_step_independent_of_end_time = true;

    /**
     * Boolean indicating if adaptive time-stepping is enabled.
     * To enable it, enable either SimulationControl::adapt_with_cfl or
     * SimulationControl::adapt_with_capillary_time_step_ratio.
     *
     * @remark By default, this is set to @p false since both
     * SimulationControl::adapt_with_cfl and
     * SimulationControl::adapt_with_capillary_time_step_ratio are set to
     * @p false by default.
     */
    bool time_step_adaptation_required = false;

    /**
     * Boolean indicating if the CFL condition should be controlling the
     * simulation time step.
     *
     * @remark By default, this is set to @p false.
     */
    bool adapt_with_cfl = false;

    // Max CFL
    double maxCFL = 1.;

    // Max time step
    double max_dt = 1e6;

    /**
     * Boolean indicating if the capillary time-step ratio should be controlling
     * the simulation time step
     *
     * @remark By default, this is set to @p false.
     */
    bool adapt_with_capillary_time_step_ratio = false;

    /**
     * The capillary time-step ratio (CTR) corresponds to the imposed value for
     * the ratio between the current time step and the capillary time-step
     * constraint (Δt/Δt_σ) to be maximally respected when
     * SimulationControl::adapt_with_capillary_time_step_ratio
     * is set to @p true.
     *
     * @remark By default, this is set to 1.
     */
    double max_capillary_time_step_ratio = 1.0;

    // Aimed tolerance at which simulation is stopped
    double stop_tolerance = 1e-10;

    // Rate of increase of the time step value
    double adaptative_time_step_scaling = 1.1;

    // BDF startup time scaling
    double startup_timestep_scaling = 0.4;

    // True if the time step should be overridden upon restart
    bool override_time_step_on_restart = false;

    // Number of mesh adaptation (steady simulations)
    unsigned int number_mesh_adaptation = 0;

    // Folder for simulation output
    std::string output_folder = "./";

    // Prefix for simulation output
    std::string output_name = "out";

    /**
     * Criterion used to end a transient simulation: either a maximum time
     * (SimulationControl::time_end) or a maximum number of transient
     * iterations (SimulationControl::iteration_end).
     *
     * @remark By default, this is set to @p EndControl::time. The default value
     * matters since this structure is also filled programmatically (e.g. in
     * unit tests) instead of being parsed from a parameter file.
     */
    enum class EndControl : std::uint8_t
    {
      iteration,
      time
    } end_control = EndControl::time;

    enum class OutputControl : std::uint8_t
    {
      iteration,
      time
    } output_control = OutputControl::iteration;

    enum class BDFStartupMethods : std::uint8_t
    {
      initial_solution,
      multiple_step_bdf,
    } bdf_startup_method = BDFStartupMethods::multiple_step_bdf;

    // Frequency of the output (for iteration output control)
    unsigned int output_iteration_frequency = 1;

    // Time frequency of the output (for time output control)
    double output_time_frequency = -1;

    // Output at specific times (for time output control)
    std::vector<double> output_times_vector{-1.};

    // Time window for file output (for both iteration and time output control)
    std::vector<double> output_time_interval{
      0.,
      std::numeric_limits<double>::max()};

    // Enable output of the boundaries
    bool output_boundaries = false;

    // Frequency of the log output to the terminal
    unsigned int log_frequency = 1;

    // Display precision of the log output to the terminal
    unsigned int log_precision = 6;

    // Subdivisions of the results in the output
    unsigned int subdivision = 1;

    // Subdivisions of the results in the output
    unsigned int group_files = 1;

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
    double T_solidus = 0.;

    // Liquidus temperature - Units in K
    double T_liquidus = 1.;

    // Latent enthalpy for the phase change - Units in J/kg
    double latent_enthalpy = 1.;

    // Specific heat of liquid - Units in J/(kg*K)
    double cp_l = 1.;

    // Specific heat of solid - Units in J/(kg*K)
    double cp_s = 1.;

    // Thermal conductivity of liquid - Units in W/(m*K)
    double thermal_conductivity_l = 1.;

    // Thermal conductivity of solid - Units in W/(m*K)
    double thermal_conductivity_s = 1.;

    // Thermal expansion coefficient of liquid - Units in 1/K
    double thermal_expansion_l = 1.;

    // Thermal expansion coefficient of solid - Units in 1/K
    double thermal_expansion_s = 0.;

    // kinematic viscosity of liquid - Units in m^2/(s)
    double kinematic_viscosity_l = 1.;

    // kinematic viscosity of solid - Units in m^2/(s)
    double kinematic_viscosity_s = 1.;

    // Darcy penalty of liquid - Units in 1/(s)
    double penalty_l = 0.;

    // Darcy penalty of solid - Units in 1/(s)
    double penalty_s = 0.;

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
    double K = 1.0;
    // Flow behavior index
    double n = 0.5;
    // Minimal shear rate magnitude for which we calculate kinematic viscosity,
    // since power-law does not allow for minimal kinematic viscosity
    double shear_rate_min = 0.001;

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
    double kinematic_viscosity_0 = 1.0;
    // Hypothetical kinematic viscosity of the flow when the shear rate is very
    // high
    double kinematic_viscosity_inf = 1.0;
    // Relaxation time
    double lambda = 1.0;
    // Carreau parameter
    double a = 2.0;
    // Power parameter
    double n = 0.5;

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
    double tracer_diffusivity_inside        = 1.0;
    double tracer_diffusivity_outside       = 1.0;
    double tracer_reaction_constant_inside  = 0.0;
    double tracer_reaction_constant_outside = 0.0;
    double thickness                        = 1.0;

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
    double tracer_diffusivity_interface       = 0.0;
    double tracer_diffusivity_bulk            = 1.0;
    double tracer_reaction_constant_interface = 0.0;
    double tracer_reaction_constant_bulk      = 0.0;
    double thickness                          = 1.0;

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
    double density_ref = 1.2;
    // Specific gas constant in J/kg/K
    double R = 287.05;
    // Absolute temperature of the ideal gas in K
    double T = 293.15;

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
    double surface_tension_coefficient = 0.0;
    // Temperature of the reference state corresponding to the surface tension
    // coefficient (T_0) in K
    double T_0 = 0.0;
    // Surface tension gradient with respect to the temperature (dsigma/dT) in
    // N/(m*K)
    double surface_tension_gradient = 0.0;

    // Solidus temperature - Units in K
    double T_solidus = 0.;

    // Liquidus temperature - Units in K
    double T_liquidus = 1.;

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
    double mobility_cahn_hilliard_constant = 1e-7;

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
    double kinematic_viscosity = 1.;
    // volumetric mass density (rho) in units of kg/m^3
    double density = 1.;
    // specific heat capacity (cp) in J/K/kg
    double specific_heat = 1.;
    // thermal conductivity (k) in W/m/K
    double thermal_conductivity = 1.;
    // thermal expansion coefficient (alpha) in 1/K
    double thermal_expansion = 1.;

    // tracer diffusivity in L^2/s
    double tracer_diffusivity = 0.;

    // tracer reaction constant in 1/s^[order]
    double tracer_reaction_constant = 0.;
    // tracer reaction order
    double tracer_reaction_order = 1.;
    // When the reaction order <1, concentrations are used in the denominator of
    // a few terms. This threshold prevents the introduction of NaN entries as
    // long as it is set above 0.
    double tracer_reaction_threshold = 1e-8;

    // Phase change parameters
    PhaseChange phase_change_parameters;

    // non-Newtonian model parameters
    enum class RheologicalModel : std::int8_t
    {
      powerlaw,
      carreau,
      newtonian,
      phase_change
    } rheological_model = RheologicalModel::newtonian;
    NonNewtonian non_newtonian_parameters;

    enum class DensityModel : std::int8_t
    {
      constant,
      isothermal_ideal_gas
    } density_model = DensityModel::constant;
    IsothermalIdealGasDensityParameters isothermal_ideal_gas_density_parameters;

    enum class SpecificHeatModel : std::int8_t
    {
      constant,
      phase_change
    } specific_heat_model = SpecificHeatModel::constant;

    enum class ThermalConductivityModel : std::int8_t
    {
      constant,
      linear,
      phase_change
    } thermal_conductivity_model = ThermalConductivityModel::constant;

    enum class ThermalExpansionModel : std::int8_t
    {
      constant,
      phase_change
    } thermal_expansion_model = ThermalExpansionModel::constant;

    enum class TracerDiffusivityModel : std::int8_t
    {
      constant,
      immersed_boundary_tanh,
      immersed_boundary_gaussian
    } tracer_diffusivity_model = TracerDiffusivityModel::constant;

    enum class TracerReactionPrefactorModel : std::int8_t
    {
      none,
      constant,
      immersed_boundary_tanh,
      immersed_boundary_gaussian
    } tracer_reaction_prefactor_model = TracerReactionPrefactorModel::constant;

    enum class ElectricConductivityModel : std::int8_t
    {
      constant,
      polynomial
    } electric_conductivity_model = ElectricConductivityModel::constant;

    enum class ElectricPermittivityModel : std::int8_t
    {
      constant,
      polynomial
    } electric_permittivity_model = ElectricPermittivityModel::constant;

    enum class MagneticPermeabilityModel : std::int8_t
    {
      constant,
      polynomial
    } magnetic_permeability_model = MagneticPermeabilityModel::constant;

    // Struct that contains the parameters to handle physical properties when
    // immersed solids tanh models are used
    ImmersedSolidTanhParameters immersed_solid_tanh_parameters;

    // Struct that contains the parameters to handle physical properties when
    // immersed solids Gaussian models are used
    ImmersedSolidGaussianParameters immersed_solid_gaussian_parameters;

    // Linear thermal conductivity parameters: k = k_A0 + k_A1 * T
    double k_A0 = 0.;
    double k_A1 = 0.;

    /// Material parameters for the electromagnetics solver. Only the
    /// permittivity and the permability are complex-valued because they can
    /// store electromagnetic energy in the medium and also dissipate it.

    /// Electric_conductivity
    double              electric_conductivity = 0.;
    std::vector<double> electric_conductivity_polynomial_coefficients;

    /// Magnetic_permeability
    double              magnetic_permeability_real = 1.;
    double              magnetic_permeability_imag = 0.;
    std::vector<double> magnetic_permeability_real_polynomial_coefficients;
    std::vector<double> magnetic_permeability_imag_polynomial_coefficients;

    /// Electric_permittivity
    double              electric_permittivity_real = 1.;
    double              electric_permittivity_imag = 0.;
    std::vector<double> electric_permittivity_real_polynomial_coefficients;
    std::vector<double> electric_permittivity_imag_polynomial_coefficients;
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
    } material_interaction_type = MaterialInteractionsType::fluid_fluid;

    // Surface tension models
    enum class SurfaceTensionModel : std::int8_t
    {
      constant,
      linear,
      phase_change
    } surface_tension_model = SurfaceTensionModel::constant;
    SurfaceTensionParameters surface_tension_parameters;

    // Cahn-Hilliard mobility models
    enum class MobilityCahnHilliardModel : std::int8_t
    {
      constant,
      quartic
    } mobility_cahn_hilliard_model = MobilityCahnHilliardModel::constant;
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
    unsigned int                  number_of_fluids = 1;
    static constexpr unsigned int max_fluids       = 2;

    // Solid objects for conjugated simulations
    std::vector<Material>         solids;
    unsigned int                  number_of_solids = 0;
    static constexpr unsigned int max_solids       = 1;

    // Fluid-fluid or fluid-solid interactions
    std::vector<MaterialInteractions> material_interactions;
    unsigned int                      number_of_material_interactions = 0;
    static constexpr unsigned int     max_material_interactions       = 3;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_fluid_interactions_with_material_interaction_ids;
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>
      fluid_solid_interactions_with_material_interaction_ids;

    /*
     * Reference Temperature for all physical properties of fluids and solids.
     * Currently, this is only used by the thermal expansion models.
     */
    double reference_temperature = 0.;

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
    bool enable = false;

    /// Total number of constraints (maximum of 1 per fluid)
    unsigned int number_of_constraints = 0;

    /// Identifiers of fluids that are constrained
    std::vector<unsigned int> fluid_ids;

    /// Absolute tolerance applied on filtered phase indicator for
    /// constrained cell selection
    std::vector<double> filtered_phase_indicator_tolerance;

    /// Lower threshold values of the constraining field (temperature)
    std::vector<double> temperature_min_values;

    /// Upper threshold values of the constraining field (temperature)
    std::vector<double> temperature_max_values;

    /// Enable/disable (@p true/false) the definition of a plane for geometrical
    /// restrictions on the domain where the stasis constraint is applied.
    bool enable_domain_restriction_with_plane = false;

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
    bool use_default_stabilization = true;

    bool heat_transfer_dcdd_stabilization = false;

    /// Boolean indicating if the DCDD stabilization term for the CLS phase
    /// fraction should be assembled (@p true) or not (@p false).
    bool cls_dcdd_stabilization = true;

    // Diffusion factor scaling the DCDD stabilization term in the CLS equation
    double dcdd_diffusion_coeff = 0.5;

    // Pressure scaling factor used to facilitate the linear solving when
    // pressure and velocity have very different scales
    double pressure_scaling_factor = 1.;

    enum class NavierStokesStabilization : std::int8_t
    {
      pspg_supg,
      gls,
      grad_div
    } stabilization = NavierStokesStabilization::pspg_supg;

    enum class ScalarLimiters : std::int8_t
    {
      moe,
      none
    } scalar_limiter = ScalarLimiters::none;

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

    Type type = Type::none;

    bool write_time_in_error_table = false;

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
    Verbosity verbosity = Verbosity::quiet;

    // Enable force post-processing
    bool calculate_force = false;

    // Enable torque post-processing
    bool calculate_torque = false;

    // Frequency of the output
    unsigned int calculation_frequency = 1;

    // Frequency of the output
    unsigned int output_frequency = 1;

    // Output precision
    unsigned int output_precision = 10;

    // Prefix for simulation output
    std::string force_output_name = "force";

    // Prefix for the torque output
    std::string torque_output_name = "torque";

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
    bool enable_radiation = false;

    // Parameters for the radiation term at the melt pool free surface
    double Stefan_Boltzmann_constant = 5.6703e-8;
    double emissivity                = 0.6;
    double Tinf                      = 0.0;

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
    bool activate_laser = false;

    // Type of laser model used in simulations. With "exponential_decay", the
    // laser acts as a volumetric source, whereas, with
    // "heat_flux_cls_interface", the laser behaves as a surface flux at the
    // interface between fluids (CLS auxiliary physic must be enabled to use
    // this model).
    enum class LaserType : std::int8_t
    {
      exponential_decay,
      gaussian_heat_flux_cls_interface,
      uniform_heat_flux_cls_interface
    } laser_type = LaserType::gaussian_heat_flux_cls_interface;

    /** Activates the angle of incidence dependence. The heat flux
     * term is multiplied by a factor \f$ \boldsymbol{n_\Gamma} \cdot
     * \boldsymbol{d_\mathrm{laser}}\f$ where \f$ \boldsymbol{n_\Gamma}\f$ is
     * the surface unit normal vector and \f$\boldsymbol{d_\mathrm{laser}}\f$ is
     * the unit direction vector of the laser.*/
    bool enable_angle_of_incidence_dependence = true;

    // Laser concentration factor indicates the definition of the beam radius.
    // In almost all the articles, it is assumed equal to 2.0
    double concentration_factor = 2.0;

    // Laser power in W
    double laser_power = 0.0;

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
    double laser_absorptivity = 0.5;

    // Penetration depth of the laser
    double penetration_depth = 1.0;

    // Laser beam radius on the melt pool surface
    double beam_radius = 0.0;

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
    } beam_orientation; // Not given a default: its default token ("z-") is
                        // only valid for dim == 3, so no single value works
                        // for both instantiations of this templated class
                        // (see rotation_axis below for the same issue).

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
    double rotation_angle = 0.0;

    // rotation axis: not given a default here, its declare_entry default
    // "0.0, 0.0, 1.0" is dim == 3 specific and only used in 3D (see
    // beam_orientation above for the same templated-dim issue)
    Tensor<1, dim> rotation_axis;

    // rotation matrix
    Tensor<2, dim> rotation_matrix;

    // Laser scan path indicates the path of the laser focal point during a
    // simulation
    std::shared_ptr<Functions::ParsedFunction<dim>> laser_scan_path;

    // Start and end time of the laser operation
    double start_time = 0.0;
    double end_time   = 1.0;

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
   *
   * @tparam dim Denotes the number of spatial dimensions.
   */
  template <int dim>
  struct PostProcessing
  {
    /**
     * @brief Isocontour bounding box parameters. This is used to get the
     * evolution of the bounding values (In 2D: \f$x_\mathrm{min}\f$,
     * \f$x_\mathrm{max}\f$, \f$y_\mathrm{min}\f$, and \f$y_\mathrm{max}\f$;
     * In 3D: \f$x_\mathrm{min}\f$, \f$x_\mathrm{max}\f$, \f$y_\mathrm{min}\f$,
     * \f$y_\mathrm{max}\f$, \f$z_\mathrm{min}\f$, and \f$z_\mathrm{max}\f$) of
     * one or many specified isocontours.
     * The data of each isocontour is written in a separate file.
     *
     * @remark At the moment, this is only implemented for the variables
     * 'temperature' and 'phase' (CLS phase indicator).
     */
    struct IsocontourBoundingBoxes
    {
      /**
       * Set of parameters related to an isocontour required for identifying the
       * bounding values (\f$x_\mathrm{min}\f$, \f$x_\mathrm{max}\f$,
       * \f$y_\mathrm{min}\f$, \f$y_\mathrm{max}\f$, \f$z_\mathrm{min}\f$, and
       * \f$z_\mathrm{max}\f$) and writing them in an output file.
       */
      struct Isocontour
      {
        /// Isocontour value
        double isovalue;

        /// Isocontour output filename
        std::string output_name;
      };

      /// Number of monitored isocontours
      unsigned int number_of_isocontour_bounding_boxes = 0;

      /// Map that regroups all isocontours of a same variable under the same
      /// variable key
      std::multimap<Variable, std::pair<unsigned int, Isocontour>>
        ids_and_isocontours_per_variable;

      /**
       * @brief Declare the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parse the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);
    };

    /**
     * @brief Probing points parameters. This is used to monitor the evolution
     * of variables (e.g., temperature) at specified points in the domain.
     *
     * The probing points are grouped according the evaluated variables in the
     * map ProbingPoints::probing_points_per_variable for easier access when
     * postprocessing in the different physics. The names (output filenames) of
     * the probes are stored in the vector
     * ProbingPoints::probing_points_output_names. The indices of the vector
     * correspond to the probe IDs.
     *
     * @remark The points defined must be part of the domain.
     */
    struct ProbingPoints
    {
      /**
       * @brief IDs and locations associated with probe points of a given
       * variable.
       */
      struct ProbingPointsPerVariable
      {
        /** Identifiers of the probe points. This is required since multiple
         * variable can be evaluated a same probing location, but all fields of
         * a same probe are written in the same file. */
        std::vector<unsigned int> ids;

        /// Point locations where we want to evaluate values
        std::vector<Point<dim>> points;
      };

      /// Number of probing points
      unsigned int number_of_probing_points = 0;

      /// Maximum number of probing points
      unsigned int max_number_of_probing_points = 25;

      /// Probing point output names ordered by ID
      std::vector<std::string> probing_points_output_names;

      /**
       * Map that regroups all probing points information (ID and coordinates)
       * of a same variable under the same variable key. This map is used within
       * the physics to quickly identify the probes that request an evaluation
       * of one of its variables.
       *
       * @remark This is the only place where the probing points are stored.
       * They are directly parsed in this map from the parameter file.
       */
      std::map<Variable, ProbingPointsPerVariable> probing_points_per_variable;

      /**
       * @brief Declares the parameters in the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      declare_parameters(ParameterHandler &prm);

      /**
       * @brief Parses the parameters from the parameter handler.
       *
       * @param[in,out] prm The parameter handler.
       */
      void
      parse_parameters(ParameterHandler &prm);

      /**
       * @brief Adds an entry to probing_points_per_variable.
       *
       * @param[in] variable Variable of the monitored quantity.
       * @param[in] id Identifier of the probe.
       * @param[in] point Location of the probe.
       */
      void inline add_probing_point(const Variable     &variable,
                                    const unsigned int &id,
                                    const Point<dim>   &point)
      {
        // Initialize new variable entry if there is none and get probes of the
        // variable
        ProbingPointsPerVariable &probes_of_variable =
          probing_points_per_variable.try_emplace(variable).first->second;

        auto &ids = probes_of_variable.ids;

        AssertThrow(std::find(ids.begin(), ids.end(), id) == ids.end(),
                    ExcMessage(
                      "For the probing point " + Utilities::int_to_string(id) +
                      ", you are specifying multiple times the variable '" +
                      get_variable_string(variable) +
                      "'. Please only specify it once."));

        // Emplace back id and point coordinates if entry is not duplicated
        ids.emplace_back(id);
        probes_of_variable.points.emplace_back(point);
      }
    };

    /// Contains all isocontour bounding boxes
    IsocontourBoundingBoxes isocontour_bounding_boxes;

    /// Contains all probing points
    ProbingPoints probing_points;

    /// Verbosity level of the post-processed quantities
    Verbosity verbosity = Verbosity::quiet;

    /// Enable total kinetic energy post-processing
    bool calculate_kinetic_energy = false;

    /// Enable total enstrophy post-processing
    bool calculate_enstrophy = false;

    /// Enable pressure power post-processing
    bool calculate_pressure_power = false;

    /// Enable viscous dissipation post-processing
    bool calculate_viscous_dissipation = false;

    /// Enable calculating apparent viscosity
    bool calculate_apparent_viscosity = false;

    /// Enable velocity post-processing
    bool calculate_average_velocities = false;

    /// Enable average temperature and average heat flux post-processing
    bool calculate_average_temp_and_hf = false;

    /// Enable pressure drop post-processing
    bool calculate_pressure_drop = false;

    /// The inlet boundary ID for pressure drop calculation
    unsigned int inlet_boundary_id = 0;

    /// The outlet boundary ID for pressure drop calculation
    unsigned int outlet_boundary_id = 1;

    /// Enable flow rate post-processing
    bool calculate_flow_rate = false;

    /// Enable tracer flow rate post-processing
    bool calculate_tracer_flow_rate = false;

    /// Set initial time to start calculations for velocities
    double initial_time_for_average_velocities = 0.0;

    /// Set initial time to start calculations for average temperature and
    /// average heat flux
    double initial_time_for_average_temp_and_hf = 0.0;

    /// Frequency of the calculation of post-processed quantities
    unsigned int calculation_frequency = 1;

    /// Frequency of the output
    unsigned int output_frequency = 1;

    /// Prefix for kinetic energy output
    std::string kinetic_energy_output_name = "kinetic_energy";

    /// Prefix for pressure drop output
    std::string pressure_drop_output_name = "pressure_drop";

    /// Prefix for flow rate output
    std::string flow_rate_output_name = "flow_rate";

    /// Prefix for tracer flow rate output
    std::string tracer_flow_rate_output_name = "tracer_flow_rate";

    /// Prefix for the enstrophy output
    std::string enstrophy_output_name = "enstrophy";

    /// Prefix for the pressure power output
    std::string pressure_power_output_name = "pressure_power";

    /// Prefix for the viscous dissipation output
    std::string viscous_dissipation_output_name = "viscous_dissipation";

    /// Prefix for the apparent viscosity output
    std::string apparent_viscosity_output_name = "apparent_viscosity";

    /// Enable tracer statistics
    bool calculate_tracer_statistics = false;

    /// Prefix for the tracer output
    std::string tracer_output_name = "tracer_statistics";

    /// Enable temperature statistics
    bool calculate_temperature_statistics = false;

    /// Prefix for the temperature output
    std::string temperature_output_name = "temperature_statistics";

    /// Enable calculation of algebraic melt volume in phase change problems
    bool calculate_algebraic_melt_volume = false;

    /// Prefix for the algebraic melt volume output
    std::string algebraic_melt_volume_output_name = "melt_volume_alge";

    /// Enable calculation of geometric melt volume in phase change problems
    bool calculate_geometric_melt_volume = false;

    /// Prefix for the geometric melt volume output
    std::string geometric_melt_volume_output_name = "melt_volume_geo";

    /// FluidIndicator corresponding to the fluid which has phase change
    Parameters::FluidIndicator monitored_fluid_with_phase_change =
      Parameters::FluidIndicator::fluid0;

    /// Melting temperature iso-value
    double melting_temperature = 0;

    /// Enable heat flux calculation
    bool calculate_heat_flux = false;

    /// Prefix for the total heat flux output
    std::string heat_flux_output_name = "heat_flux";

    /// Fluid domain, used when post-processing a multiphase simulation
    Parameters::FluidIndicator postprocessed_fluid =
      Parameters::FluidIndicator::both;

    /// Enable barycenter calculation for fluid 1 in CLS and Cahn-Hilliard
    /// simulations
    bool calculate_barycenter = false;

    /// Prefix for the CLS and Cahn-Hilliard barycenter output
    std::string barycenter_output_name = "barycenter_information";

    /// Enable phase statistics
    bool calculate_phase_statistics = false;

    /// Prefix for the phase output
    std::string phase_output_name = "phase_statistics";

    /// Enable mass conservation calculation for both fluids in CLS simulations
    bool calculate_mass_conservation = true;

    /// Prefix for the CLS mass conservation output
    std::string mass_conservation_output_name = "mass_conservation_information";

    /// Enable energies calculation on the domain in Cahn-Hilliard simulations
    bool calculate_phase_energy = false;

    /// Prefix for the energy output in Cahn-Hilliard simulations
    std::string phase_energy_output_name = "phase_energy";

    /// Enable calculation of total fluid volume and total particles volume in
    /// cfd-dem simulation
    bool calculate_phase_volumes = false;

    /// prefix for the total volume output in cfd-dem simulation
    std::string phase_volumes_output_name = "phase_volumes";

    /// Enable output of Q-criterion field
    bool output_q_criterion = true;

    /// Enable output of vorticity field
    bool output_vorticity = true;

    /// Enable output of velocity gradient field
    bool output_velocity_gradient = true;

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief FEM - The finite element section
   * controls the properties of the finite element method. This section
   * controls the interpolation degree of polynomial integration
   * and the number of quadrature points within the cells.
   */
  struct FEM
  {
    // Interpolation degree velocity
    unsigned int velocity_degree = 1;

    // Interpolation degree pressure
    unsigned int pressure_degree = 1;

    // Interpolation degree void fraction
    unsigned int void_fraction_degree = 1;

    // Interpolation degree temperature
    unsigned int temperature_degree = 1;

    // Interpolation degree tracer
    unsigned int tracer_degree = 1;

    // Switch tracer to DG formulation instead of CG
    bool tracer_uses_dg = false;

    // Interpolation degree cls model
    unsigned int CLS_degree = 1;

    // Switch cls to DG formulation instead of CG
    bool CLS_uses_dg = false;

    // Interpolation degree Cahn-Hilliard
    unsigned int phase_cahn_hilliard_degree     = 1;
    unsigned int potential_cahn_hilliard_degree = 1;

    // Option for bubble enrichment functions
    bool enable_bubble_function_velocity = false;
    bool enable_bubble_function_pressure = false;

    /// Polynomial degree for the different electromagnetics spaces.
    /// The trial degree sets the polynomial degree for the solution and the
    /// test degree is used for the computation of the cell matrix necessary
    /// when solving a system using the DPG method.
    unsigned int electromagnetics_trial_degree = 1;
    unsigned int electromagnetics_test_degree  = 2;

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

    Verbosity verbosity = Verbosity::verbose;

    // Type of non-linear solver
    SolverType solver = SolverType::newton;

    // Kinsol solver strategy
    KinsolStrategy kinsol_strategy = KinsolStrategy::line_search;

    // Tolerance
    double tolerance = 1e-6;

    // Maximal number of iterations for the Newton solver
    unsigned int max_iterations = 10;

    // Force RHS recalculation at the beginning of every non-linear steps
    // This is required if there is a fixed point component to the non-linear
    // solver that is changed at the beginning of every newton iteration.
    // This is notably the case of the sharp edge method.
    // The default value of this parameter is false.
    bool force_rhs_calculation = false;

    // Matrix reconstruction tolerance
    // This parameter controls the reconstruction of the system matrix
    // If the residual after a newton step is lower than previous_residual *
    // matrix_tolerance then that iteration is considered sufficient and the
    // matrix is not reassembled at the next iteration.
    double matrix_tolerance = 0.1;

    // Relative Tolerance
    double step_tolerance = 0.9;

    // Carry jacobian matrix over to the new non-linear problem
    bool reuse_matrix = false;

    // Reuse preconditioner for the next non-linear iterations
    bool reuse_preconditioner = false;

    // Abort solver if non-linear solution has not reached tolerance
    bool abort_at_convergence_failure = false;

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

    SolverType solver = SolverType::gmres;

    /// Verbosity of linear solver
    Verbosity verbosity = Verbosity::verbose;

    /// Flag to rescale linear and non-linear residuals by the sqrt of the
    /// triangulation volume
    bool rescale_residual_by_volume = false;

    /// Relative residuals of the iterative solver
    double relative_residual = 1e-3;

    /// Minimum residual of the iterative solver
    double minimum_residual = 1e-12;

    /// Maximum number of iterations
    int max_iterations = 1000;

    /// Maximum number of krylov vectors
    int max_krylov_vectors = 100;

    /// Enable hessians in jacobian
    bool enable_hessians_jacobian = true;

    /// Enable hessians in residual
    bool enable_hessians_residual = true;

    /// Type of preconditioner
    enum class PreconditionerType : std::int8_t
    {
      ilu,  /// Incomplete LU Factorization
      amg,  /// Algebraic Multigrid
      lsmg, /// Local-Smoothing Geometric Multigrid
      gcmg, /// Global-Coarsening Multigrid
      none  /// No preconditioner
    };
    PreconditionerType preconditioner = PreconditionerType::ilu;

    /// ILU fill
    unsigned int ilu_precond_fill = 0;

    /// ILU absolute tolerance
    double ilu_precond_atol = 1e-12;

    /// ILU relative tolerance
    double ilu_precond_rtol = 1.00;

    /// AMG parameters either as linear solver preconditioner or as
    /// preconditioner of a coarse-grid solver for LSMG or GCMG

    /// ILU fill for smoother
    unsigned int amg_precond_ilu_fill = 0;

    /// ILU absolute tolerance for smoother
    double amg_precond_ilu_atol = 1e-12;

    /// ILU relative tolerance for smoother
    double amg_precond_ilu_rtol = 1.00;

    /// AMG aggregation threshold
    double amg_aggregation_threshold = 1e-14;

    /// AMG number of cycles
    unsigned int amg_n_cycles = 1;

    /// AMG W_cycle
    bool amg_w_cycles = false;

    /// AMG Smoother sweeps
    unsigned int amg_smoother_sweeps = 2;

    /// AMG Smoother overlap
    unsigned int amg_smoother_overlap = 1;

    /// Block linear solver to throw error.
    bool force_linear_solver_continuation = false;

    /// MG min level
    int mg_min_level = -1;

    /// MG minimum number of cells per level
    int mg_level_min_cells = -1;

    /// MG intermediate level
    int mg_int_level = -1;

    /// MG enable hessians in jacobian
    bool mg_enable_hessians_jacobian = true;

    /// Type of multigrid
    enum class MultigridCoarseningSequenceType : std::int8_t
    {
      h,
      p,
      hp,
      ph
    };
    MultigridCoarseningSequenceType mg_coarsening_type =
      MultigridCoarseningSequenceType::h;

    /// Type of p coarsening sequence
    MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
      mg_p_coarsening_type = MGTransferGlobalCoarseningTools::
        PolynomialCoarseningSequenceType::decrease_by_one;

    /// Minimum polynomial degree for p coarsening sequence
    unsigned int mg_p_min_coarsening_degree = 1;

    /// MG smoother number of iterations
    int mg_smoother_iterations = 10;

    /// MG smoother relaxation parameter
    double mg_smoother_relaxation = 0.5;

    /// Type of preconditioner for the MG smoother
    enum class MultigridSmootherPreconditionerType : std::int8_t
    {
      InverseDiagonal,
      AdditiveSchwarzMethod,
      Chebyshev
    };
    MultigridSmootherPreconditionerType mg_smoother_preconditioner_type =
      MultigridSmootherPreconditionerType::InverseDiagonal;

    /// Chebyshev smoother polynomial degree
    unsigned int mg_smoother_chebyshev_degree = 3;

    /// Chebyshev smoother smoothing range (lambda_max / lambda_min)
    double mg_smoother_chebyshev_smoothing_range = 15;

    /// Chebyshev smoother eigenvalue-estimation CG/Lanczos iterations
    unsigned int mg_smoother_chebyshev_eig_cg_n_iterations = 10;

    /// MG eigenvalue estimation for smoother relaxation parameter
    bool mg_smoother_eig_estimation = true;

    /// MG smoothing range to set range between eigenvalues
    double eig_estimation_smoothing_range = 10;

    /// MG number of cg iterations to find eigenvalue
    int eig_estimation_cg_n_iterations = 10;

    /// MG print max, min, eigenvalues
    Verbosity eig_estimation_verbose = Verbosity::verbose;

    /// Type of coarse grid solver
    enum class CoarseGridSolverType : std::int8_t
    {
      gmres,
      amg,
      ilu,
      direct
    };
    CoarseGridSolverType mg_coarse_grid_solver = CoarseGridSolverType::direct;

    /// MG use FE_Q_iso_Q1 elements for coarse grid
    bool mg_use_fe_q_iso_q1 = false;

    /// MG coarse-grid solver maximum number of iterations
    int mg_gmres_max_iterations = 2000;

    /// MG coarse-grid solver tolerance
    double mg_gmres_tolerance = 1e-14;

    /// MG coarse-grid solver reduce
    double mg_gmres_reduce = 1e-4;

    /// MG coarse-grid solver maximum number of krylov vectors
    int mg_gmres_max_krylov_vectors = 30;

    /// MG coarse-grid solver preconditioner
    PreconditionerType mg_gmres_preconditioner = PreconditionerType::amg;

    /// MG use default parameters for AMG
    bool mg_amg_use_default_parameters = false;

    /// MG information about levels
    Verbosity mg_verbosity = Verbosity::verbose;

    static void
    declare_parameters(ParameterHandler &prm, const std::string &physics_name);
    void
    parse_parameters(ParameterHandler &prm, const std::string &physics_name);
  };

  /**
   * @brief Mesh - Parameters that control mesh reading and mesh generation.
   */
  template <int dim, int spacedim = dim>
  class Mesh
  {
  public:
    // GMSH or dealii
    enum class Type : std::int8_t
    {
      gmsh,
      dealii,
      lethe
    };
    Type type = Type::dealii;

    /// File name of the mesh
    std::string file_name = "none";

    /// Name of the grid in GridTools
    std::string grid_type = "hyper_cube";

    /// Arguments of the GridTools
    std::string grid_arguments = "-1 : 1 : false";

    /// Initial refinement level of primitive mesh
    unsigned int initial_refinement = 0;

    /// Initial refinement level of primitive mesh near user-defined boundary
    /// conditions
    unsigned int initial_refinement_at_boundaries = 0;

    /// List of boundary ids to refine
    std::vector<int> boundaries_to_refine;

    /// Enable fixing initial refinement from a target size
    bool refine_until_target_size = false;

    /// Allow the use of a simplex mesh
    bool simplex = false;

    /// Target size when automatically refining initial mesh
    double target_size = 1.;

    /// Enable checking the input grid for diamond-shaped cells
    bool check_for_diamond_cells = false;

    /* A boolean parameter which enables adding the neighbor boundary cells of
    * boundary cells in DEM simulations. This parameter should only be enabled
    * for simulations with concave geometries (for instance particles inside a
    * drum). In simulations with convex geometries, it must not be enabled.
    * This is also reported to users in a warning in
     find_boundary_cells_information.*/
    bool expand_particle_wall_contact_search = false;

    // Grid displacement at initiation
    Tensor<1, spacedim> translation;

    // Grid rotation at initiation
    Tensor<1, spacedim> rotation_axis = [] {
      Tensor<1, spacedim> axis;
      axis[0] = 1.;
      return axis;
    }();
    double rotation_angle = 0.;

    /// Rescale the grid by the scale factor
    double scale = 1.;

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
    /// Error estimator for the variable
    enum class ErrorEstimator : std::int8_t
    {
      kelly,
      dpg
    } error_estimator = ErrorEstimator::kelly;

    // Coarsening fraction
    double coarsening_fraction = 0.05;

    // Refinement fraction
    double refinement_fraction = 0.1;
  };

  /**
   * @brief MeshAdaption - Parameters that control dynamic mesh adaptation.
   * Dynamic mesh adaptation in Lethe is very flexible and can be both local
   * and global.
   */
  struct MeshAdaptation
  {
    // Initial adaptive refinement
    unsigned int initial_refinement = 0;

    // Type of mesh adaptation
    enum class Type : std::int8_t
    {
      none,
      uniform,
      adaptive
    } type = Type::none;

    // Map containing the refinement variables
    std::map<Variable, MultipleAdaptationParameters> variables;
    // declaration for parsing variables
    Variable                     vars = Variable::velocity;
    MultipleAdaptationParameters var_adaptation_param;

    /// Decision factor for adaptive refinement (number or fraction)
    enum class FractionType : std::int8_t
    {
      number,
      fraction
    } fractionType = FractionType::number;

    // Maximum number of elements
    unsigned int maximum_number_elements = 100000000;

    // Maximum refinement level
    unsigned int maximum_refinement_level = 10;

    // Minimum refinement level
    unsigned int minimum_refinement_level = 0;

    // Refinement after frequency iter
    unsigned int frequency               = 1;
    bool         refinement_at_frequency = true;

    // Enable the control of the mesh refinement to target a specific number of
    // elements equal to the maximum number of elements.
    bool mesh_controller_is_enabled = false;

    // Specifies if mesh adaptation can be used on certain boundaries
    bool is_boundary_refinement_fixed = false;

    // List of boundary ids to fix to their initial refinement state
    std::vector<int> boundaries_to_fix;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * Container of the parameters for box refinements. The regions to refine can
   * be described by either a GMSH or a deal.II mesh.
   */
  template <int dim, int spacedim = dim>
  struct MeshBoxRefinement
  {
    /// Number of boxes delimiting refinement regions
    unsigned int number_of_refinement_boxes = 0;

    /**
     * Maximum number of refinement boxes that can be defined by a user.
     *
     * @remark A maximal value has to be initialized to fill the vectors with
     * parameter declarations.
     */
    const unsigned int max_number_of_refinement_boxes = 20;

    /**
     * Shared pointer of a vector of GMSH and deal.II meshes representing
     * refinement areas.
     */
    std::shared_ptr<std::vector<Mesh<dim, spacedim>>> refinement_boxes_meshes =
      std::make_shared<std::vector<Mesh<dim, spacedim>>>(
        max_number_of_refinement_boxes);

    /// Vector of additional refinement values of the different boxes
    std::vector<unsigned int> box_additional_refinements =
      std::vector<unsigned int>(max_number_of_refinement_boxes);

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
    bool enabled = false;

    enum class TestType : std::int8_t
    {
      particles,
      mobility_status,
      subdomain
    } test_type = TestType::particles;

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
    std::string  filename   = "restart";
    bool         restart    = false;
    bool         checkpoint = false;
    unsigned int frequency  = 1;
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
    enum class PermeabilityModel : std::int8_t
    {
      none,               ///< No Darcy source term
      darcy_phase_change, /**< Phase change Darcy source term which applies a
                    linear penalization depending on the phase change model */
      carman_kozeny_phase_change /**< Phase change Carman-Kozeny source term
                    which applies a non-linear penalization depending on the
                    phase change model */
    };

    enum class RotatingFrameType : std::int8_t
    {
      none,
      srf
    };

    RotatingFrameType rotating_frame_type = RotatingFrameType::none;
    double            omega_x             = 0.;
    double            omega_y             = 0.;
    double            omega_z             = 0.;

    /// Type of permeability model applied to the Navier-Stokes equations
    PermeabilityModel permeability_model = PermeabilityModel::none;

    /// Indicates which fluids are with liquid-solid phase change.
    FluidIndicator fluid_with_phase_change = FluidIndicator::fluid0;

    /**
     * Enable the multiplication of the Darcy force term
     * (\f$\vec{F}_\mathrm{Darcy}\f$) by the density for dimensional consistency
     * when solving the pressure rather than the kinematic pressure in the
     * momentum balance.
     *
     * \f[
     * \vec{F}_\mathrm{Darcy} = \rho K \vec{u}
     * \f]
     *
     * with \f$\rho\f$ the density, \f$K\f$ the Darcy penalty, and \f$\vec{u}\f$
     * the velocity.
     */
    bool enable_darcy_multiply_by_density = false;

    /// Permeability area of the pseudo-porous bed (solid phase).
    std::vector<double> carman_kozeny_permeability_area;

    /// Tolerance in the Carman-Kozeny source term that avoids division by zero.
    std::vector<double> carman_kozeny_tolerance;

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
    unsigned int nb_particles = 1;
    // Boolean to determine whether the Navier-Stokes equations are
    // solved inside the particles.
    bool assemble_navier_stokes_inside = false;

    // Polynomial degree of the IB stencil
    unsigned int stencil_degree = 2;
    // The length ratio used for the stencil calculation of the IB condition.
    double length_ratio = 4;
    // Boolean controlling whether extrapolation is used to impose the
    // immersed boundary condition. If false, the IB condition is directly
    // imposed using nearest neighbors. the immersed boundary condition or not.
    // If it is set to false, all cut cells are fully imposed on the IB.
    bool enable_extrapolation = true;

    // Boolean for the calculation of the force at the IB
    bool calculate_force_ib = true;
    // Boolean for extra vtu field output
    bool enable_extra_sharp_interface_vtu_output_field = false;
    // The name of the output file for the forces on the IB.
    std::string ib_force_output_file = "ib_force";
    // Particles pvd file name
    std::string ib_particles_pvd_file = "ib_particles_data";
    // Boolean for printing DEM information
    bool print_dem = true;

    // Number of initial refinements around each particle
    unsigned int initial_refinement = 0;
    // Boolean to enable distance-based coarsening away from particles
    bool enable_coarsening = false;
    // Inner distance factor of the refinement zone
    double refinement_inside_distance_factor = 0.5;
    // Outer distance factor of the refinement zone
    double refinement_outside_distance_factor = 1.5;
    // Outer distance factor beyond which distance-based coarsening is
    // allowed
    double coarsening_distance_factor = 2.0;
    // Boolean for time-dependent mesh refinement according to the particle's
    // current position the refinement zone.
    bool time_extrapolation_of_refinement_zone = false;

    // Number of DEM time steps per CFD time step.
    unsigned int coupling_frequency = 1000;
    // Relaxation parameter for the CFD-DEM coupling.
    double alpha = 1;
    // Frequency at which the contact search is performed at the CFD time
    // scale (once every X CFD time steps)
    int contact_search_frequency = 1;
    // Particles' radius multiplier used to calculate the effective radius of
    // contact search
    double contact_search_radius_factor = 3;
    // Boolean for lubrication force
    bool enable_lubrication_force = true;
    // The maximal range for which the lubrication force is evaluated. This
    // variable multiplies the smallest cell diameter to obtain the actual
    // range.
    double lubrication_range_max = 2;
    // The minimal range for which the lubrication force is evaluated. This
    // variable multiplies the smallest cell diameter to obtain the actual
    // range.
    double lubrication_range_min = 0.1;
    // Tolerance for the particle dynamics nonlinear solver.
    double particle_nonlinear_tolerance = 1e-6;
    // Boolean for the explicit calculation of the contact impulse in the
    // CFD-DEM coupling impulsion in the CFD-DEM coupling.
    bool explicit_contact_impulsion_calculation = false;
    // Boolean for explicit evaluation of the particle's position in the CFD-DEM
    // coupling
    bool explicit_position_integration_calculation = false;
    // Boolean for approximation of the contact radius. If true, the effective
    // radius replaces the actual particle's local curvature radius in the
    // calculation.
    bool approximate_radius_for_contact = false;
    // Function defining the gravitational acceleration vector used by the
    // CFD-DEM calculation. Given a default so that declare_parameters() can
    // safely dereference it on a default-constructed instance.
    std::shared_ptr<Functions::ParsedFunction<dim>> f_gravity =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    // Young's modulus of the wall
    double wall_youngs_modulus = 100000000;
    // Poisson ratio of the wall
    double wall_poisson_ratio = 0.3;
    // Rolling friction coefficient of the wall
    double wall_rolling_friction_coefficient = 0;
    // Sliding friction coefficient of the wall
    double wall_friction_coefficient = 0;
    // Coefficient of restitution of the wall
    double wall_restitution_coefficient = 1;

    // Boolean for loading particles from an independent file. If true, the
    // definition of the particles in the Particle subsection in the .prm file
    // is ignored. of the particle subsection.
    bool load_particles_from_file = false;
    // Name of the independent file containing particles' information at
    // insertion. Only used if load_particles_from_file is true.
    std::string particles_file = "particles";
  };

  /**
   * @brief FlowControl - Set average velocity on a boundary (CFD) or the domain
   * (CFD-DEM).
   */
  struct DynamicFlowControl
  {
    // Enable flow control
    bool enable_flow_control = false;

    // Average velocity target (L/t)
    double average_velocity_0 = 0.;

    // Boundary id at flow inlet
    unsigned int boundary_flow_id = 0;

    // Flow direction (x=0, y=1 ,z=2)
    unsigned int flow_direction = 0;

    // Initial beta
    double beta_0 = 0.;

    // Relaxation coefficient for beta force controller
    // beta_n+1 = beta_n + alpha * (...)
    double alpha = 1.;

    // If beta at n+1 step is in this threshold over beta at n step, beta n+1
    // is kept as beta n. This avoids a new term of force in the matrix
    double beta_threshold = 0.0;

    // Type of verbosity for the flow control
    Verbosity verbosity = Verbosity::quiet;

    // Apply scaled beta force to particles
    bool enable_beta_particle = false;

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
    } evaporative_mass_flux_model_type = EvaporativeMassFluxModelType::constant;

    bool enable_evaporation_cooling = false;
    bool enable_recoil_pressure     = false;

    // Parameters for the evaporation terms at the melt pool free surface
    double evaporation_mass_flux       = 0.0;
    double evaporation_coefficient     = 0.82;
    double recoil_pressure_coefficient = 0.56;
    double molar_mass                  = 1.0;
    double boiling_temperature         = 1.0;
    double latent_heat_evaporation     = 0.0;
    double ambient_pressure            = 101325;
    double ambient_gas_density         = 1.0;
    double liquid_density              = 10.0;
    double universal_gas_constant      = 8.3145;

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
    bool enable = false;
    /// Type of mortar interface
    enum class InterfaceType : std::int8_t
    {
      circular,
      linear
    } interface_type = InterfaceType::circular;
    /// Mesh parameters for the rotor part
    std::shared_ptr<Mesh<dim>> rotor_mesh;
    /// Boundary ID # of the rotor at the rotor-stator interface
    unsigned int rotor_boundary_id = 1;
    /// Boundary ID # of the stator at the rotor-stator interface
    unsigned int stator_boundary_id = 2;
    /// Center of rotation of the rotor domain
    Point<dim> center_of_rotation;
    /// Rotation axis direction
    unsigned int rotation_axis_direction = 2;
    /// Rotation angle of the rotor domain in radians
    std::shared_ptr<Functions::ParsedFunction<dim>> rotor_rotation_angle;
    /// Angular velocity of the rotor domain
    std::shared_ptr<Functions::ParsedFunction<dim>> rotor_angular_velocity;
    /// Penalty factor for mortar elements
    double sip_factor = 1.;
    /// Oversampling factor for quadrature points
    unsigned int oversampling_factor = 1;
    /// Tolerance used for rotor-stator interface radius computation
    double radius_tolerance = 1e-8;
    /// Cell weight for load balancing of cells with mortar interfaces
    unsigned int cell_weight = 1000;
    /// Type of verbosity for mortar
    Verbosity verbosity = Verbosity::quiet;

    void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
  };


} // namespace Parameters
#endif
