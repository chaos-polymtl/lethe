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
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#include <solvers/initial_conditions.h>

namespace Parameters
{
  void
  Ramp_n::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ramp n");
    {
      prm.declare_entry(
        "n0",
        "1.0",
        Patterns::Double(),
        "First n value with which to start the initial condition");

      prm.declare_entry(
        "n_iter",
        "5",
        Patterns::Double(),
        "Number of iterations uses in the ramp before reaching the final n value");

      prm.declare_entry("alpha",
                        "0.5",
                        Patterns::Double(),
                        "Coefficient used for n-spacing.");
    }
    prm.leave_subsection();
  }


  void
  Ramp_n::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ramp n");
    {
      n0     = prm.get_double("n0");
      n_iter = prm.get_double("n_iter");
      alpha  = prm.get_double("alpha");
    }
    prm.leave_subsection();
  }

  extern template class InitialConditions<2>;
  extern template class InitialConditions<3>;
} // namespace Parameters
