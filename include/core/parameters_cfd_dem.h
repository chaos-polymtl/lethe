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

#include <deal.II/base/parsed_function.h>

#include <core/parameters.h>


using namespace dealii;
/**
 * The analytical solution class provides an interface for all common
 * dlement required for the calculation of analytical solution
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
    dem
  };

  enum class DragModel
  {
    difelice,
    rong
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
    std::string                    dem_file_name;
    double                         l2_smoothing_factor;
    double                         l2_lower_bound;
    double                         l2_upper_bound;
    unsigned int                   void_fraction_fem_degree;
  };

  struct CFDDEM
  {
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    bool         shock_capturing;
    bool         grad_div;
    bool         full_stress_tensor;
    double       reference_velocity;
    DragModel    drag_model;
    bool         post_processing;
    unsigned int inlet_boundary_id;
    unsigned int outlet_boundary_id;
  };



} // namespace Parameters

#endif
