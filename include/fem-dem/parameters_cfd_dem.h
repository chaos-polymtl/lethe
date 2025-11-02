// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_parameters_cfd_dem_h
#define lethe_parameters_cfd_dem_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;
/**
 * The analytical solution class provides an interface for all common
 * elementS required for the calculation of analytical solution
 * All equation-specific analytical solution should derive
 * from the base class but also call it's declare_parameters and
 *parse_parameters routine. This allows specialize class to focus on their
 *specificity and forget about other non-specific elements that are generic to
 *the calculation of analytical solutions
 **/
namespace Parameters
{
  enum class VoidFractionMode
  {
    function,
    pcm, // The particle centered method
    qcm, // The quadratured centered method
    spm  // The satellite point method (divided approach)
  };

  enum class VoidFractionQuadratureRule
  {
    gauss,        // Apply gauss quadrature rule (default)
    gauss_lobatto // Apply gauss-lobatto quadrature rule
  };

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
   *   The drag is computed the previous velocity of the particles and the
   * fluid. The resulting drag is applied to the particles. On the fluid, the
   * drag is used to calculate the momentum transfer coefficient and the new
   * fluid velocity velocity is calculated using an implicit formulation.
   * Provides extensive stability, but for larger time-step this may lead to a
   * slight violation of Newton's third law since the drag applied on the fluid
   *   is not strictly equal to the drag applied to the particles.
   *
   * - **fully_explicit**:
   *   The drag is evaluated entirely using known quantities from the previous
   *   time step. The outcome conditionally stable and may require very small
   *   coupling time-step when the density of the fluid is much lower than the
   * density of the particles
   *
   */
  enum class DragCoupling
  {
    fully_implicit,
    semi_implicit,
    fully_explicit
  };


  enum class VANSModel
  {
    modelA,
    modelB
  };


  template <int dim>
  class VoidFractionParameters
  {
  public:
    VoidFractionParameters()
      : void_fraction(1)
    {}

    /**
     * Default destructor.
     */
    virtual ~VoidFractionParameters() = default;

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);

  public:
    VoidFractionMode               mode;
    Functions::ParsedFunction<dim> void_fraction;
    bool                           read_dem;
    std::string                    dem_file_name;
    double                         l2_smoothing_length;
    unsigned int                   particle_refinement_factor;
    double                         qcm_sphere_diameter;
    bool                           qcm_sphere_equal_cell_volume;
    VoidFractionQuadratureRule     quadrature_rule;
    unsigned int                   n_quadrature_points;
    bool                           project_particle_velocity;
  };

  struct CFDDEM
  {
    bool         grad_div;
    DragModel    drag_model;
    DragCoupling drag_coupling;
    VANSModel    vans_model;
    unsigned int coupling_frequency;
    bool         drag_force;
    bool         buoyancy_force;
    bool         shear_force;
    bool         pressure_force;
    bool         saffman_lift_force;
    bool         magnus_lift_force;
    bool         rotational_viscous_torque;
    bool         vortical_viscous_torque;
    bool         void_fraction_time_derivative;
    bool         interpolated_void_fraction;
    double       cstar;
    bool         implicit_stabilization;
    bool         particle_statistics;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters
#endif
