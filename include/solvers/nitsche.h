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
 * Author: Carole-Anne Daunais, Polytechnique Montreal, 2020 -
 */


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <core/parameters.h>

#ifndef nitsche_h
#define nitsche_h

using namespace dealii;

namespace Parameters
{
  template <int dim>
  class Nitsche
  {
  public:

    Nitsche()
     // : solid_velocity(std::make_shared<Functions::ParsedFunction<dim>>(dim))
      {}

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
    
    // Solid mesh
    Parameters::Mesh solid_mesh;

    // Penalization term
    double beta;

    // Solid velocity
    std::shared_ptr<Functions::ParsedFunction<dim>> solid_velocity;
  };
  
  template <int dim>
  void
  Nitsche<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("nitsche");
    {
      solid_mesh.declare_parameters(prm);
      prm.declare_entry("beta",
                        "1",
                        Patterns::Double(),
                        "Penalization term for Nitsche method");
      prm.enter_subsection("solid velocity");
      solid_velocity->declare_parameters(prm, dim);
      if (dim == 2)
        prm.set("Function expression", "0; 0");
      if (dim == 3)
        prm.set("Function expression", "0; 0; 0");
      prm.leave_subsection();

    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Nitsche<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("nitsche");
    {
      solid_mesh.parse_parameters(prm);
      beta = prm.get_double("beta");
      prm.enter_subsection("solid velocity");
        solid_velocity->parse_parameters(prm);
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  
} // namespace Parameters

#endif
