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
 * Author: Toni EL Geitani, Bruno Blais, Polytechnique Montreal, 2020-
 */

#ifndef lethe_parameters_cfd_dem_h
#define lethe_parameters_cfd_dem_h

#include <core/parameters.h>
#include <core/parameters_lagrangian.h>

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

  enum class DragModel
  {
    difelice,
    rong,
    dallavalle,
    kochhill,
    beetstra,
    gidaspow
  };

  enum class VANSModel
  {
    modelA,
    modelB
  };


  template <int dim>
  class VoidFraction
  {
  public:
    VoidFraction()
      : void_fraction(1)

    {}

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);


  public:
    VoidFractionMode               mode;
    Functions::ParsedFunction<dim> void_fraction;
    bool                           read_dem;
    bool                           bound_void_fraction;
    std::string                    dem_file_name;
    double                         l2_smoothing_factor;
    double                         l2_lower_bound;
    double                         l2_upper_bound;
    unsigned int                   particle_refinement_factor;
    double                         qcm_sphere_diameter;
    bool                           qcm_sphere_equal_cell_volume;
  };

  struct CFDDEM
  {
    bool         grad_div;
    DragModel    drag_model;
    VANSModel    vans_model;
    bool         post_processing;
    unsigned int inlet_boundary_id;
    unsigned int outlet_boundary_id;
    unsigned int coupling_frequency;
    bool         drag_force;
    bool         buoyancy_force;
    bool         shear_force;
    bool         pressure_force;
    bool         saffman_lift_force;
    bool         void_fraction_time_derivative;
    bool         interpolated_void_fraction;
    double       cstar;
    bool         implicit_stabilization;

    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters

#endif
