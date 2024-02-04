/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 */


#ifndef lethe_tracer_drift_velocity_h
#define lethe_tracer_drift_velocity_h

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/lac/vector.h>

#include <memory>

using namespace dealii;
/**
 * Drift velocity that is added to the tracer velocity to account for simplified multiphase simulations
 * using drift-flux modeling.
 **/

namespace Parameters
{
  template <int dim>
  class TracerDriftVelocity
  {
  public:
    TracerDriftVelocity()
    {
      drift_velocity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
    }

    virtual void
    declare_parameters(ParameterHandler &prm);
    virtual void
    parse_parameters(ParameterHandler &prm);

    // Drift velocity
    std::shared_ptr<Functions::ParsedFunction<dim>> drift_velocity;
  };

  template <int dim>
  void
  TracerDriftVelocity<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("tracer drift velocity");

    prm.enter_subsection("drift velocity");
    drift_velocity->declare_parameters(prm, dim);
    prm.leave_subsection();

    prm.leave_subsection();
  }

  template <int dim>
  void
  TracerDriftVelocity<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("tracer drift velocity");

    prm.enter_subsection("drift velocity");
    drift_velocity->parse_parameters(prm);
    prm.leave_subsection();

    prm.leave_subsection();
  }
}


#endif
