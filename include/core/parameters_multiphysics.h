// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/*
 * This file defines the parameters in the parameter namespace
 * that pertain to multiphysics simulations
 */

#ifndef lethe_parameters_multiphysics_h
#define lethe_parameters_multiphysics_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief Different interface reinitialization method types:
   *  - none
   *  - projection_based_sharpening: projection-based interface sharpening
   *  - pde_based: PDE-based reinitialization
   *  - geometric: geometric redistanciation
   */
  enum class ReinitializationMethodType : std::int8_t
  {
    none,
    projection_based_sharpening,
    pde_based,
    geometric
  };

  /**
   * @brief Different projection-based interface sharpening types:
   *  - constant: the sharpening threshold is the same throughout the
   * simulation,
   *  - adaptive: the sharpening threshold is determined by binary search, to
   * ensure mass conservation of the monitored phase
   */
  enum class SharpeningType : std::int8_t
  {
    constant,
    adaptive
  };

  /**
   * @brief Different transformation function types for the signed distance:
   *  - tanh: hyperbolic tangent
   *  - piecewise_polynomial: 4th degree piecewise polynomial. The latter takes
   *    the following form:
   *
   *    if \f$d < 0 \f$
   *    \f$\phi = 0.5 - 0.5(4d' + 6d'^2 + 4*d'^3 + d'^4)\f$
   *
   *    else
   *    \f$\phi = 0.5 - 0.5(4d' - 6d'^2 + 4*d'^3 - d'^4)\f$
   *
   *    with \f$d' = d/d_\mathrm{max}\f$, where \f$\phi\f$ is the phase
   *    fraction, \f$d\f$ the signed distance and \f$d_\mathrm{max}\f$ the
   *    maximum redistanciation distance. This transformation clamps the phase
   *    fraction to 0 or 1 when \f$d = \pm d_\mathrm{max}\f$.
   */
  enum class RedistanciationTransformationType : std::int8_t
  {
    tanh,
    piecewise_polynomial
  };

  /**
   * @brief Different phase indicator filtering types:
   * - none: no filter wil be applied on the calculated phase indicator
   * - tanh: the tanh filter function will be applied to the phase indicator,
   * a \f$\beta\f$ parameter influencing the interface definition must be
   * defined
   */
  enum class FilterType : std::int8_t
  {
    none,
    clip,
    tanh
  };

  enum class EpsilonSetMethod : std::int8_t
  {
    automatic,
    manual
  };

  /**
   * @brief Verbosity options for the epsilon parameter.
   */
  enum class EpsilonVerbosity : std::int8_t
  {
    /// Epsilon related information will not be displayed on terminal
    quiet,
    /// Epsilon value will be displayed on terminal for every steady and
    /// transient iteration.
    verbose
  };

  /**
   * @brief Types of waveguide port modes:
   * - Transverse Electric (TE): the electric field is transverse to the
   * direction of propagation.
   * - Transverse Magnetic (TM): the magnetic field is transverse to the
   * direction of propagation.
   */
  enum class WaveguideMode : std::int8_t
  {
    /// Transverse Electric mode
    TE,
    /// Transverse Magnetic mode
    TM
  };

  /**
   * @brief Type of electromagnetic scaling to apply to the solution of the time-harmonic Maxwell solver after solving the linear system. This is relevant when the user wants to recover the physical solution in dimensional units from  the dimensionless solution.
   * - none: no scaling will be applied to the solution after solving the linear
   * system, so the solution will be in dimensionless units.
   * - electric_field: the solution will be scaled by the electric scaling
   * factor given from the user input parameter. This factor needs to match the
   * dimensionality of the problem.
   * - magnetic_field: the solution will be scaled by the magnetic scaling
   * factor given from the user input parameter. This factor needs to match the
   * dimensionality of the problem.
   * - power: the solution will be scaled by the power given from the user input
   * parameter in waveguide inlets. This option is only relevant when there is
   * at least one waveguide inlet in the problem. Note that it is the user's
   * responsibility to ensure that if the problem also applies non-zero electric
   * or magnetic field Dirichlet boundary conditions, those need to be scaled
   * accordingly.
   */
  enum class ElectromagneticScalingType : std::int8_t
  {
    none,
    electric_field,
    magnetic_field,
    power
  };

  /**
   * @brief Parameters for the time coupling strategy of the time-harmonic Maxwell solver with the other physics solvers.
   * - none: the time-harmonic Maxwell solver will only be solved at the first
   * iteration of the time integration and will not be solved in the subsequent
   * iterations.
   * - iteration: the time-harmonic Maxwell solver will be solved every N time
   * iterations, where N is given by the user input parameter "coupling
   * iteration".
   * - time: the time-harmonic Maxwell solver will be solved when crossing a
   * time multiple, e.g. every T seconds, where the time multiple is given by
   * the user input parameter "coupling time".
   * - threshold: the time-harmonic Maxwell solver will be solved when the
   * maximum change in the material properties of the problem change beyond a
   * certain threshold, where the threshold value is given by the user input
   * parameter "coupling threshold".
   */
  enum class TimeHarmonicMaxwellCouplingStrategy : std::int8_t
  {
    none,
    iteration,
    time,
    threshold
  };

  /**
   * @brief CahnHilliard_PhaseFilter - Defines the parameters for the phase filtration of CahnHilliard physics
   */
  struct CahnHilliard_PhaseFilter
  {
    // Type of filter
    Parameters::FilterType type = FilterType::none;

    // \f$beta\f$ value for the tanh filter
    double beta = 20;

    // Type of verbosity for the phase filter
    Parameters::Verbosity verbosity = Parameters::Verbosity::quiet;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief CLS_InterfaceSharpening - Defines the parameters for
   * interface sharpening in the CLS solver.
   */
  struct CLS_InterfaceSharpening
  {
    // Interface sharpening parameters. The sharpening method and parameters are
    // explained in the dam break CLS example:
    // https://chaos-polymtl.github.io/lethe/examples/multiphysics/dam-break-CLS/dam-break-CLS.html

    bool enable = false;

    Parameters::SharpeningType type = SharpeningType::constant;

    // Parameters for constant sharpening
    double threshold = 0.5;

    // Parameters for adaptive sharpening
    double threshold_max_deviation = 0.20;
    int    max_iterations          = 20;

    // Other sharpening parameters
    double interface_sharpness = 2;

    bool monitoring = false;

    /// Conservation tolerance on the fluid monitored,
    /// used with adaptive sharpening
    double tolerance = 1e-6;

    Parameters::FluidIndicator monitored_fluid =
      Parameters::FluidIndicator::fluid_1;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for the calculation of surface tension
   * force in the CLS solver.
   */
  struct CLS_SurfaceTensionForce
  {
    bool enable = false;

    double phase_indicator_gradient_diffusion_factor = 4;
    double curvature_diffusion_factor                = 1;

    bool output_cls_auxiliary_fields = false;

    // Type of verbosity for the surface tension force calculation
    Parameters::Verbosity verbosity = Parameters::Verbosity::quiet;

    // Enable marangoni effect
    bool enable_marangoni_effect = false;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for the phase filtration
   */
  struct CLS_PhaseFilter
  {
    // Type of filter
    Parameters::FilterType type = FilterType::none;

    // $$\beta$$ value for the tanh filter
    double beta = 20;

    // Type of verbosity for the phase filter
    Parameters::Verbosity verbosity = Parameters::Verbosity::quiet;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for PDE-based reinitialization of the interface
   * used with the CLS solver.
   */
  struct CLS_PDEBasedInterfaceReinitialization
  {
    /// Enables/Disables the PDE-based interface reinitialization.
    bool enable = false;
    /**
     * Enables/Disables @p pvtu format outputs of the PDE-based interface
     * reinitialization steps of the last simulated time step.
     * The files are stored in a folder named
     * @p pde-based-reinitialization-steps-output located inside the
     * <tt>output path</tt> folder specified in the <tt>simulation control</tt>
     * subsection.
     * */
    bool output_reinitialization_steps = false;
    /// Constant multiplying the mesh-size in the evaluation of the diffusion
    /// coefficient.
    double diffusivity_multiplier = 1.;
    /// Constant representing the power to which the mesh-size is elevated in
    /// the evaluation of the diffusion coefficient.
    double diffusivity_power = 1.;
    /// Constant factor used in the computation of the artificial time step of
    /// the reinitialization equation
    double dtau_factor = 0.25;
    /// Steady-state criterion used for the artificial time-stepping scheme.
    double steady_state_criterion = 1e-4;
    /// Maximum number of reinitialization steps.
    double max_steps_number = 10000;

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
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for geometric reinitialization of the interface
   * used with the CLS solver.
   */
  struct CLS_GeometricInterfaceReinitialization
  {
    /// Enables/Disables the geometric interface reinitialization.
    bool enable = false;
    /// Enables/Disables the output of the signed distance field
    bool output_signed_distance = false;
    /// Maximum reinitialization distance value
    double max_reinitialization_distance = 1.;
    /// Transformation type transforming the signed distance to a phase
    /// indicator
    RedistanciationTransformationType transformation_type =
      RedistanciationTransformationType::tanh;
    /// Interface thickness for the tanh transformation
    double tanh_thickness = 1.;

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
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for interface reinitialization methods
   * used within the CLS solver. It stores the parameters for the three
   * available methods (projection-, PDE-, and geometric based
   * reinitialization).
   */
  struct CLS_ReinitializationMethod
  {
    /// Reinitialization method type
    Parameters::ReinitializationMethodType reinitialization_method_type =
      ReinitializationMethodType::none;

    /// Reinitialization frequency at every \f$x\f$ time steps the CLS phase
    /// fraction field will be regularized
    int frequency = 10;

    /// Type of verbosity of the PDE-based interface reinitialization solver.
    Parameters::Verbosity verbosity = Parameters::Verbosity::quiet;

    /// Interface sharpening parameters
    Parameters::CLS_InterfaceSharpening sharpening;

    /// PDE-based interface reinitialization parameters
    Parameters::CLS_PDEBasedInterfaceReinitialization
      pde_based_interface_reinitialization;

    /// Geometric interface reinitialization parameters
    Parameters::CLS_GeometricInterfaceReinitialization
      geometric_interface_reinitialization;

    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    void
    declare_parameters(ParameterHandler &prm) const;

    /**
     * @brief Parse the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief CLS - Defines the parameters for free surface simulations
   * using the CLS method.
   * Has to be declared before member creation in Multiphysics structure.
   */
  struct CLS
  {
    Parameters::CLS_SurfaceTensionForce    surface_tension_force;
    Parameters::CLS_PhaseFilter            phase_filter;
    Parameters::CLS_ReinitializationMethod reinitialization_method;

    Parameters::FluidIndicator viscous_dissipative_fluid =
      Parameters::FluidIndicator::fluid_1;

    // artificial diffusivity (diffusion coefficient) (in L^2/s) added to the
    // CLS transport equation. This parameter is zero by default, and can be
    // increased to improve the wetting of the phases in the vicinity of
    // boundaries
    double diffusivity = 0;

    bool compressible = false;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct CahnHilliard
  {
    // Smoothing parameter \f$\xi\f$ in the Cahn-Hilliard equations (potential
    // equation)
    double potential_smoothing_coefficient = 1;

    // Epsilon set strategy (automatic|manual)
    Parameters::EpsilonSetMethod epsilon_set_method =
      EpsilonSetMethod::automatic;

    // Epsilon verbosity
    Parameters::EpsilonVerbosity epsilon_verbosity = EpsilonVerbosity::quiet;

    // Epsilon value in the Cahn-Hilliard equations
    double epsilon = 1.0;

    // Phase filtration parameters
    Parameters::CahnHilliard_PhaseFilter cahn_hilliard_phase_filter;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Parameters of the phase-space evolution of the QMOM moments, which
   * is delegated to QMOMKit.
   */
  struct QMOMPhaseEvolution
  {
    /// Enable the adaptive number of quadrature nodes
    bool enable_adaptive_node_count = false;

    /**
     * @brief Declare the "phase evolution" subsection.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm) const;

    /**
     * @brief Parse the "phase evolution" subsection.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Operator splitting scheme between the physical-space transport (T)
   * and the phase-space evolution (P) of the QMOM fields:
   *  - lie: T then P over the whole time step (first order);
   *  - strang: P over half a time step, T, then P over half a time step
   *    (second order).
   */
  enum class QMOMSplittingScheme : std::int8_t
  {
    lie,
    strang
  };

  /**
   * @brief Process parameters of the lipid nanoparticle (LNP) formulation
   * simulated with the QMOM physics: operating conditions, properties of the
   * aqueous and organic phases and shape factor of the particles. All of them
   * are in SI units.
   */
  struct QMOMPharma
  {
    /* Operating conditions */
    /// Operating temperature [K]
    double operating_temperature = 295.;
    /// Ratio of the aqueous phase flow rate to the organic phase flow rate [-]
    double flow_rate_ratio = 1.;
    /// Total flow rate of the aqueous and organic phases [m^3/s]
    double total_flow_rate = 1.;
    /// Mixing time scale of the aqueous and organic phases [s]
    double mixing_time_scale = 0.;
    /// Initial lipid concentration in the organic phase [mol/m^3]
    double initial_concentration_organic = 1.;
    /// Ion concentration in the mixture [mol/m^3]
    double ion_concentration_mixture = 1.;

    /* Aqueous phase */
    /// Density of the aqueous phase [kg/m^3]
    double density_aqueous = 997.;
    /// Molar mass of the aqueous phase [kg/mol]
    double molar_mass_aqueous = 0.018015;
    /// Mole fraction solubility of the lipid in the pure aqueous phase [-]
    double mole_fraction_solubility_aqueous = 2e-8;
    /// Relative permittivity of the aqueous phase [-]
    double permittivity_aqueous = 78.2;

    /* Organic phase */
    /// Density of the organic phase [kg/m^3]
    double density_organic = 789.;
    /// Molar mass of the organic phase [kg/mol]
    double molar_mass_organic = 0.046069;
    /// Mole fraction solubility of the lipid in the pure organic phase [-]
    double mole_fraction_solubility_organic = 2.9e-3;
    /// Relative permittivity of the organic phase [-]
    double permittivity_organic = 24.3;

    /* Mixture and particles */
    /// Dynamic viscosity of the aqueous and organic mixture [Pa s]
    double viscosity_mixture = 2e-3;
    /// Volume shape factor of the particles, pi/6 for spheres [-]
    double shape_factor = std::numbers::pi / 6.;

    /**
     * @brief Declare the "pharma" subsection.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm) const;

    /**
     * @brief Parse the "pharma" subsection.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Physical constants used by the QMOM physics, in SI units. Their
   * default values are the exact values of the 2019 redefinition of the SI
   * base units and the CODATA 2022 value of the vacuum permittivity.
   */
  struct QMOMPhysicalConstants
  {
    /// Boltzmann constant [J/K]
    double boltzmann_constant = 1.380649e-23;
    /// Avogadro number [1/mol]
    double avogadro_number = 6.02214076e23;
    /// Elementary charge [C]
    double elementary_charge = 1.602176634e-19;
    /// Vacuum permittivity [F/m]
    double permittivity_vacuum = 8.8541878188e-12;

    /**
     * @brief Declare the "physical constants" subsection.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm) const;

    /**
     * @brief Parse the "physical constants" subsection.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Family of the particle size distribution from which the moments of
   * the QMOM physics are computed. The raw moments of each family are computed
   * analytically by QMOMKit:
   *  - none: no particles, all the moments are zero;
   *  - normal: mean and standard deviation;
   *  - lognormal: mean and standard deviation of the logarithm of the size;
   *  - gamma: shape and scale;
   *  - weibull: shape and scale;
   *  - beta: shape parameters alpha and beta, on the interval [0, 1];
   *  - chi: number of degrees of freedom;
   *  - chi_squared: number of degrees of freedom.
   */
  enum class QMOMDistributionType : std::int8_t
  {
    none,
    normal,
    lognormal,
    gamma,
    weibull,
    beta,
    chi,
    chi_squared
  };

  /**
   * @brief Number of parameters of a particle size distribution family, not
   * counting the number density.
   *
   * @param[in] distribution_type Family of the distribution.
   *
   * @return Number of parameters of the family.
   */
  unsigned int
  get_number_of_distribution_parameters(
    const QMOMDistributionType distribution_type);

  /**
   * @brief Functions of space and time prescribing the fields of the QMOM
   * physics. They are used for the initial conditions, the Dirichlet boundary
   * conditions, the analytical solution and the source terms.
   *
   * The scalars are given as a list of expressions, one per scalar. The
   * moments are not given one by one: they are the raw moments of a particle
   * size distribution, defined by its family and by a list of expressions
   * whose first entry is the number density, in particles per m^3, followed by
   * the parameters of the family.
   *
   * The number of scalars is not known when the parameters are declared. The
   * scalars function therefore accepts either a single expression, applied to
   * every scalar, or one expression per scalar. The QMOM physics verifies the
   * number of expressions once the number of scalars is known.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  struct QMOMFunctions
  {
    /// Scalars, or a single function applied to all of them
    std::shared_ptr<Functions::ParsedFunction<dim>> scalars;

    /// Family of the particle size distribution
    QMOMDistributionType distribution_type = QMOMDistributionType::none;

    /// Number density, in particles per m^3, followed by the parameters of
    /// the particle size distribution. Ignored if the distribution type is
    /// none.
    std::shared_ptr<Functions::ParsedFunction<dim>> distribution;

    /**
     * @brief Declare the scalars function in the current subsection and the
     * "distribution" subsection. All the functions are zero and the
     * distribution type is none by default.
     *
     * @param[in,out] prm Parameter handler.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the functions from the current subsection. The number of
     * components of each function is the number of expressions separated by
     * semicolons. The number of expressions of the distribution must be one
     * (the number density) plus the number of parameters of its family.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /**
     * @brief Set the time of all the functions.
     *
     * @param[in] time New time.
     */
    void
    set_time(const double time) const;
  };

  /**
   * @brief Parameters of the quadrature method of moments (QMOM) physics,
   * declared in the "QMOM" subsection.
   */
  struct QMOM
  {
    /// Number of scalars transported together with the moments
    unsigned int number_of_scalars = 2;

    /// Names of the scalars, used for the output. scalar_0, scalar_1, ... if
    /// they are not given in the parameter file.
    std::vector<std::string> scalar_names = {"scalar_0", "scalar_1"};

    /// Number of transported raw moments M_0 to M_{number_of_moments - 1}
    unsigned int number_of_moments = 2;

    /**
     * @brief Get the number of nodes of the quadrature that uses all the
     * transported moments. N nodes are defined by 2N moments, so an odd number
     * of moments requires one more node.
     *
     * @return Number of moments divided by 2, rounded up.
     */
    unsigned int
    get_number_of_nodes() const
    {
      return (number_of_moments + 1) / 2;
    }

    /// Operator splitting scheme between the physical-space transport and
    /// the phase-space evolution
    QMOMSplittingScheme splitting_scheme = QMOMSplittingScheme::lie;

    /// Include the nucleation source
    bool enable_nucleation = true;
    /// Include the coalescence source
    bool enable_coalescence = true;
    /// Include the growth source
    bool enable_growth = true;

    /// Phase-space evolution ("phase evolution" subsection)
    QMOMPhaseEvolution phase_evolution;

    /// Process parameters of the LNP formulation ("pharma" subsection)
    QMOMPharma pharma;

    /// Physical constants ("physical constants" subsection)
    QMOMPhysicalConstants physical_constants;

    /**
     * @brief Declare the "QMOM" subsection.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm) const;

    /**
     * @brief Parse the "QMOM" subsection.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };


  /**
   * @brief TimeHarmonicMaxwell - Defines the parameters for
   * time-harmonic Maxwell simulations.
   */
  template <int dim>
  struct TimeHarmonicMaxwell
  {
    // Integration method for the time coupling of the time-harmonic Maxwell
    // solver
    Parameters::TimeHarmonicMaxwellCouplingStrategy time_coupling_strategy =
      TimeHarmonicMaxwellCouplingStrategy::none;

    // Coupling parameter for the time coupling strategy. If the time coupling
    // strategy is based on the iteration count, this parameter represents the
    // number of time iterations between two consecutive resolutions of the
    // electromagnetic fields.
    unsigned int coupling_iteration = 1;

    // If the time coupling strategy is based on time, this parameter represents
    // the real time interval between two consecutive resolutions of the
    // electromagnetic fields.
    double coupling_time = 1;

    // If the time coupling strategy is based on a threshold, this parameter
    // represents the change in the electromagnetic properties of the medium
    // (i.e., permittivity, permeability or conductivity) that triggers the
    // recomputation of the electromagnetic fields.
    double coupling_threshold = 0.1;

    // We use vectors in the following to be able to define multiple waveguides
    // in the same simulation.
    unsigned int number_of_waveguide_inlets = 0;

    // Frequency of the electromagnetic wave (in Hz)
    double electromagnetic_frequency = 1;

    // Boundary ids where waveguide ports are applied so the port can be linked
    // to the right boundary
    std::vector<int> waveguide_boundary_ids;

    // Waveguide power (in W, time-averaged) to be able to link the amplitude
    // of the electromagnetic wave to the power injected in the waveguide
    std::vector<double> waveguide_power;

    // Electric field amplitude used for the normalization of the solution (in
    // V/m)
    double electric_field_amplitude = 1;

    // Magnetic field amplitude used for the normalization of the solution (in
    // A/m)
    double magnetic_field_amplitude = 1;

    // Type of scaling to apply to the solution after solving the linear system
    Parameters::ElectromagneticScalingType electromagnetic_scaling_type =
      ElectromagneticScalingType::none;

    // Waveguide mode to simulate (TE|TM)
    std::vector<Parameters::WaveguideMode> waveguide_mode;

    // Waveguide modes order (m,n) to simulate
    std::vector<unsigned int> mode_order_m;
    std::vector<unsigned int> mode_order_n;

    // Waveguide corners locations
    std::vector<std::array<Tensor<1, dim>, Utilities::fixed_power<2>(dim - 1)>>
      waveguide_corners;

    void
    declare_parameters(ParameterHandler &prm) const;

    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Multiphysics - the parameters for multiphysics simulations
   * and handles sub-physics parameters.
   */
  template <int dim>
  struct Multiphysics
  {
    bool fluid_dynamics           = true;
    bool heat_transfer            = false;
    bool tracer                   = false;
    bool CLS                      = false;
    bool cahn_hilliard            = false;
    bool electromagnetics         = false;
    bool continuum_solid_dynamics = false;
    bool qmom                     = false;

    // subparameters for heat_transfer
    bool viscous_dissipation    = false;
    bool thermal_buoyancy_force = false;
    bool microwave_heating      = false;

    Parameters::CLS                      cls_parameters;
    Parameters::CahnHilliard             cahn_hilliard_parameters;
    Parameters::TimeHarmonicMaxwell<dim> time_harmonic_maxwell_parameters;
    Parameters::QMOM                     qmom_parameters;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };
} // namespace Parameters
#endif
