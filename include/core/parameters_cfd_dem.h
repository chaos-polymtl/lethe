// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_parameters_cfd_dem_h
#define lethe_parameters_cfd_dem_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

/**
 * @file parameters_cfd_dem.h
 * @brief Parameter structures for CFD-DEM coupled simulations.
 *
 * This file defines the parameter classes and enumerations used to configure
 * CFD-DEM coupled simulations in Lethe, including void fraction computation
 * methods, drag models, drag coupling strategies, VANS model selection, and
 * the main CFDDEM parameter structure.
 */
namespace Parameters
{
  /**
   * @brief Method used to compute the void fraction field.
   */
  enum class VoidFractionMode
  {
    /// Void fraction prescribed by an analytical function.
    function,
    /// Particle centered method.
    pcm,
    /// Quadrature centered method.
    qcm,
    /// Satellite point method (divided approach).
    spm
  };

  /**
   * @brief Quadrature rule used for void fraction integration.
   */
  enum class VoidFractionQuadratureRule
  {
    /// Gauss quadrature rule (default).
    gauss,
    /// Gauss-Lobatto quadrature rule.
    gauss_lobatto
  };

  /**
   * @brief Drag force model for particle-fluid interaction in CFD-DEM.
   */
  enum class DragModel
  {
    difelice,
    rong,
    dallavalle,
    kochhill,
    beetstra,
    gidaspow
  };


  /**
   * @enum DragCoupling
   * @brief Defines the numerical coupling strategy used for drag force computation between fluid and particles.
   *
   * This enumeration specifies how the drag interaction term is treated in time
   * integration schemes for coupled CFD-DEM or multiphase simulations. The
   * choice of coupling affects both stability and computational cost. The
   * fully_implicit and fully_explicit schemes are named this way since explicit
   * is a reserved keyword in C++. Thus, we use fully_implicit and
   * fully_explicit to ensure a coherent terminology.
   *
   * - **fully_implicit**:
   *   The drag force is evaluated using the current (implicit) fluid
   *   velocities and previous (explicit) particle velocity. The resulting
   *   drag is applied to both the particles and the fluid.
   *   Provides maximum stability, but requires more Newton iterations.
   *
   * - **semi_implicit**:
   *   The drag force on the particles is computed using the previous velocities
   *   of both the particles and the fluid. When transferring momentum to the
   *   fluid, the previous fluid velocity is used to evaluate the momentum
   *   transfer coefficient, while the current fluid velocity is used to
   *   calculate relative velocity between the particle and the fluid. This
   *   provides extensive stability, but for larger time step may lead to a
   *   slight violation of Newton's third law since the drag applied on the
   *   fluid is not strictly equal to the drag applied to the particles.
   *
   * - **fully_explicit**:
   *   The drag is evaluated entirely using known quantities from the previous
   *   time step. The outcome conditionally stable and may require very small
   *   coupling time step when the density of the fluid is much lower than the
   *   density of the particles
   *
   */
  enum class DragCoupling
  {
    fully_implicit,
    semi_implicit,
    fully_explicit
  };


  /**
   * @brief Volume-Averaged Navier-Stokes (VANS) model formulation.
   */
  enum class VANSModel
  {
    /// Model A: volume-averaged continuity and momentum equations.
    modelA,
    /// Model B: alternative VANS formulation.
    modelB
  };


  /**
   * @brief Parameters controlling the void fraction computation in CFD-DEM
   * simulations.
   *
   * This class stores the method, quadrature settings, and auxiliary parameters
   * used to compute the void fraction field from the particle distribution.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class VoidFractionParameters
  {
  public:
    VoidFractionParameters()
      : void_fraction(1)
    {}

    /**
     * @brief Default destructor.
     */
    virtual ~VoidFractionParameters() = default;

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    virtual void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    virtual void
    parse_parameters(ParameterHandler &prm);

  public:
    /// Void fraction computation method.
    VoidFractionMode mode;

    /// Analytical void fraction function (used when mode is function).
    Functions::ParsedFunction<dim> void_fraction;

    /// Read particle data from a DEM simulation checkpoint.
    bool read_dem;

    /// File name of the DEM simulation checkpoint.
    std::string dem_file_name;

    /// Smoothing length for L2 projection of the void fraction.
    double l2_smoothing_length;

    /// Refinement factor applied to the particle radius for void fraction.
    unsigned int particle_refinement_factor;

    /// Sphere diameter used in the quadrature centered method.
    double qcm_sphere_diameter;

    /// Use a sphere volume equal to the cell volume in QCM.
    bool qcm_sphere_equal_cell_volume;

    /// Quadrature rule used for void fraction integration.
    VoidFractionQuadratureRule quadrature_rule;

    /// Number of quadrature points used for void fraction integration.
    unsigned int n_quadrature_points;

    /// Project particle velocity onto the fluid mesh.
    bool project_particle_velocity;
  };

  /**
   * @brief Parameters for the CFD-DEM coupled solver.
   *
   * This structure stores the coupling parameters, force model selections,
   * and stabilization options used in CFD-DEM simulations.
   */
  struct CFDDEM
  {
    /// Enable grad-div stabilization.
    bool grad_div;

    /// Drag model used for particle-fluid interaction.
    DragModel drag_model;

    /// Numerical coupling strategy for drag force computation.
    DragCoupling drag_coupling;

    /// Volume-Averaged Navier-Stokes model formulation.
    VANSModel vans_model;

    /// Frequency of DEM-CFD coupling (in DEM time steps).
    unsigned int coupling_frequency;

    /// Enable drag force on particles.
    bool drag_force;

    /// Enable buoyancy force on particles.
    bool buoyancy_force;

    /// Enable shear force on particles.
    bool shear_force;

    /// Enable pressure gradient force on particles.
    bool pressure_force;

    /// Enable Saffman lift force on particles.
    bool saffman_lift_force;

    /// Enable Magnus lift force on particles.
    bool magnus_lift_force;

    /// Enable rotational viscous torque on particles.
    bool rotational_viscous_torque;

    /// Enable vortical viscous torque on particles.
    bool vortical_viscous_torque;

    /// Include the void fraction time derivative in the equations.
    bool void_fraction_time_derivative;

    /// Use interpolated void fraction instead of cell-averaged values.
    bool interpolated_void_fraction;

    /// Stabilization constant for the void fraction equation.
    double cstar;

    /// Enable implicit stabilization of the void fraction.
    bool implicit_stabilization;

    /// Enable output of particle statistics.
    bool particle_statistics;

    /// Project particle forces onto the fluid mesh.
    bool project_particle_forces;

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters
#endif
