// SPDX-FileCopyrightText: Copyright (c) 2019, 2021-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/analytical_solutions.h>

namespace AnalyticalSolutions
{
  template <int dim>
  void
  AnalyticalSolution<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("analytical solution");
    prm.declare_entry(
      "enable",
      "false",
      Patterns::Bool(),
      "Enable the calculation of the analytical solution and L2 error");
    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether from the post-processing values should be printed "
      "Choices are <quiet|verbose>.");

    prm.declare_entry(
      "filename",
      "L2Error",
      Patterns::FileName(),
      "File name for the output for the L2Error table with respect to time or mesh ");


    {
      prm.enter_subsection("uvwp");
      uvwp.declare_parameters(prm, dim + 1);
      prm.leave_subsection();
    }
    {
      prm.enter_subsection("temperature");
      temperature.declare_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("tracer");
      tracer.declare_parameters(prm);
      prm.leave_subsection();
    }
    {
      prm.enter_subsection("VOF");
      phase.declare_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.declare_parameters(prm, 2);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("electromagnetics");
      electromagnetics.declare_parameters(prm, 4 * dim);
      prm.leave_subsection();
    }

    prm.leave_subsection();
  }

  template <int dim>
  void
  AnalyticalSolution<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("analytical solution");
    enable               = prm.get_bool("enable");
    filename             = prm.get("filename");
    const std::string op = prm.get("verbosity");
    if (op == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    if (op == "quiet")
      verbosity = Parameters::Verbosity::quiet;

    {
      prm.enter_subsection("uvwp");
      uvwp.parse_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("temperature");
      temperature.parse_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("tracer");
      tracer.parse_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("VOF");
      phase.parse_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("cahn hilliard");
      cahn_hilliard.parse_parameters(prm);
      prm.leave_subsection();
    }

    {
      prm.enter_subsection("electromagnetics");
      electromagnetics.parse_parameters(prm);
      prm.leave_subsection();
    }

    prm.leave_subsection();
  }


} // namespace AnalyticalSolutions
// Pre-compile the 2D and 3D
template class AnalyticalSolutions::AnalyticalSolution<2>;
template class AnalyticalSolutions::AnalyticalSolution<3>;
