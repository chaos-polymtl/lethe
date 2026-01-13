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

using namespace dealii;

namespace Parameters
{
  /**
   * @brief Different interface regularization method types:
   *  - none
   *  - sharpening: projection-based interface sharpening
   *  - algebraic: PDE-based reinitialization
   *  - geometric: geometric redistanciation
   */
  enum class RegularizationMethodType : std::int8_t
  {
    none,
    sharpening,
    algebraic,
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
   * @brief Different phase fraction filtering types:
   * - none: no filter wil be applied on the calculated phase fraction
   * - tanh: the tanh filter function will be applied to the phase fraction,
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
   * @brief CahnHilliard_PhaseFilter - Defines the parameters for the phase filtration of CahnHilliard physics
   */
  struct CahnHilliard_PhaseFilter
  {
    // Type of filter
    Parameters::FilterType type;

    // \f$beta\f$ value for the tanh filter
    double beta;

    // Type of verbosity for the phase filter
    Parameters::Verbosity verbosity;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief VOF_InterfaceSharpening - Defines the parameters for
   * interface sharpening in the VOF solver.
   */
  struct VOF_InterfaceSharpening
  {
    // Interface sharpening parameters. The sharpening method and parameters are
    // explained in the dam break VOF example:
    // https://chaos-polymtl.github.io/lethe/examples/multiphysics/dam-break-VOF/dam-break-VOF.html

    bool enable = false;

    Parameters::SharpeningType type;

    // Parameters for constant sharpening
    double threshold;

    // Parameters for adaptive sharpening
    double threshold_max_deviation;
    int    max_iterations;

    // Other sharpening parameters
    double interface_sharpness;

    bool monitoring;

    /// Conservation tolerance on the fluid monitored,
    /// used with adaptive sharpening
    double tolerance;

    Parameters::FluidIndicator monitored_fluid;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for the calculation of surface tension
   * force in the VOF solver.
   */
  struct VOF_SurfaceTensionForce
  {
    bool enable;

    double phase_fraction_gradient_diffusion_factor;
    double curvature_diffusion_factor;

    bool output_vof_auxiliary_fields;

    // Type of verbosity for the surface tension force calculation
    Parameters::Verbosity verbosity;

    // Enable marangoni effect
    bool enable_marangoni_effect;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for the phase filtration
   */
  struct VOF_PhaseFilter
  {
    // Type of filter
    Parameters::FilterType type;

    // $$\beta$$ value for the tanh filter
    double beta;

    // Type of verbosity for the phase filter
    Parameters::Verbosity verbosity;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Parameters for algebraic reinitialization of the interface
   * used with the VOF solver.
   */
  struct VOF_AlgebraicInterfaceReinitialization
  {
    /// Enables/Disables the algebraic interface reinitialization.
    bool enable = false;
    /**
     * Enables/Disables @p pvtu format outputs of the algebraic interface
     * reinitialization steps of the last simulated time-step.
     * The files are stored in a folder named
     * @p algebraic-reinitialization-steps-output located inside the
     * <tt>output path</tt> folder specified in the <tt>simulation control</tt>
     * subsection.
     * */
    bool output_reinitialization_steps;
    /// Constant multiplying the mesh-size in the evaluation of the diffusion
    /// coefficient.
    double diffusivity_multiplier;
    /// Constant representing the power to which the mesh-size is elevated in
    /// the evaluation of the diffusion coefficient.
    double diffusivity_power;
    /// Constant factor used is the computation of the artificial time-step of
    /// the reinitialization equation
    double dtau_factor;
    /// Steady-state criterion used for the artificial time-stepping scheme.
    double steady_state_criterion;
    /// Maximum number of reinitialization steps.
    double max_steps_number;

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
   * used with the VOF solver.
   */
  struct VOF_GeometricInterfaceReinitialization
  {
    /// Enables/Disables the geometric interface reinitialization.
    bool enable = false;
    /// Enables/Disables the output of the signed distance field
    bool output_signed_distance;
    /// Maximum reinitialization distance value
    double max_reinitialization_distance;
    /// Transformation type transforming the signed distance to a phase fraction
    RedistanciationTransformationType transformation_type;
    /// Interface thickness for the tanh transformation
    double tanh_thickness;

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
   * @brief Parameters for interface regularization methods
   * used within the VOF solver. It stores the parameters for the three
   * available methods (projection-, algebraic-, and geometric based
   * regularization).
   */
  struct VOF_RegularizationMethod
  {
    /// Regularization method type
    Parameters::RegularizationMethodType regularization_method_type;

    /// Regularization frequency at every \f$x\f$ time-steps the VOF phase
    /// fraction field will be regularized
    int frequency;

    /// Type of verbosity of the algebraic interface reinitialization solver.
    Parameters::Verbosity verbosity;

    /// Interface sharpening parameters
    Parameters::VOF_InterfaceSharpening sharpening;

    /// Algebraic interface reinitialization parameters
    Parameters::VOF_AlgebraicInterfaceReinitialization
      algebraic_interface_reinitialization;

    /// Geometric interface reinitialization parameters
    Parameters::VOF_GeometricInterfaceReinitialization
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
   * @brief VOF - Defines the parameters for free surface simulations
   * using the VOF method.
   * Has to be declared before member creation in Multiphysics structure.
   */
  struct VOF
  {
    Parameters::VOF_SurfaceTensionForce  surface_tension_force;
    Parameters::VOF_PhaseFilter          phase_filter;
    Parameters::VOF_RegularizationMethod regularization_method;

    Parameters::FluidIndicator viscous_dissipative_fluid;

    // artificial diffusivity (diffusion coefficient) (in L^2/s) added to the
    // VOF transport equation. This parameter is zero by default, and can be
    // increased to improve the wetting of the phases in the vicinity of
    // boundaries
    double diffusivity;

    bool compressible;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm);
  };

  struct CahnHilliard
  {
    // Smoothing parameter \f$\xi\f$ in the Cahn-Hilliard equations (potential
    // equation)
    double potential_smoothing_coefficient;

    // Epsilon set strategy (automatic|manual)
    Parameters::EpsilonSetMethod epsilon_set_method;

    // Epsilon verbosity
    Parameters::EpsilonVerbosity epsilon_verbosity;

    // Epsilon value in the Cahn-Hilliard equations
    double epsilon;

    // Phase filtration parameters
    Parameters::CahnHilliard_PhaseFilter cahn_hilliard_phase_filter;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };

  /**
   * @brief Multiphysics - the parameters for multiphysics simulations
   * and handles sub-physics parameters.
   */
  struct Multiphysics
  {
    bool fluid_dynamics;
    bool heat_transfer;
    bool tracer;
    bool VOF;
    bool cahn_hilliard;
    bool electromagnetics;

    // subparameters for heat_transfer
    bool viscous_dissipation;
    bool buoyancy_force;

    Parameters::VOF          vof_parameters;
    Parameters::CahnHilliard cahn_hilliard_parameters;

    void
    declare_parameters(ParameterHandler &prm) const;
    void
    parse_parameters(ParameterHandler &prm, const Dimensionality &dimensions);
  };
} // namespace Parameters
#endif
