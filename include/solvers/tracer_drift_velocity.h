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
 * Drift velocity that is added to the tracer velocity to account for simplified
 *multiphase simulations using drift-flux modeling.
 **/

namespace Parameters
{
  /**
* @brief Implements a drift velocity function to account for simplified multiphase flows.
*
* @tparam dim An integer that denotes the dimension of the space in which
* the flow is solved.

  * The drift velocity provides a simple way to model dilute disperse multiphase
flow through the tracer physics.
  **/

  template <int dim>
  class TracerDriftVelocity
  {
  public:
    TracerDriftVelocity()
    {
      drift_velocity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
    }

    /**
     * @brief Declares the parameters required by the drift velocity.
     *
     * @param prm ParameterHandler used to declare the parameters.
     *
     */
    virtual void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parses the parameters required by the analytical solution
     * within a parameter file
     *
     * @param prm ParameterHandler used to parse the parameters.
     *
     */
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
} // namespace Parameters


#endif
