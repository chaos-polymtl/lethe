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
 */
#ifndef lethe_ale_h
#define lethe_ale_h

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/lac/vector.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief An interface for all the common parameters required to consider
   * a moving frame of reference in a simplified Arbitrary Lagrangian Eulerian
   * (ALE) module.
   *
   * @tparam dim An integer that denotes the dimension of the space in which
   * the problem is solved.
   */
  template <int dim>
  class ALE
  {
  public:
    /**
     * @brief Construct a new ALE object.
     */
    ALE()
    {
      velocity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
    }


    /**
     * @brief Declare the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    virtual void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters.
     *
     * @param[in,out] prm The ParameterHandler.
     */
    virtual void
    parse_parameters(ParameterHandler &prm);


    /**
     * @brief ALE velocity function
     */
    std::shared_ptr<Functions::ParsedFunction<dim>> velocity;

    /**
     * @brief Return if the ALE module is enabled.
     *
     * @return enable A boolean that indicates if the ALE module is enabled.
     */
    bool
    enabled() const
    {
      return enable;
    }

  private:
    /**
     * @brief Boolean indicating if the ALE module is enabled.
     */
    bool enable;
  };

  template <int dim>
  void
  ALE<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ALE");
    prm.declare_entry(
      "enable",
      "false",
      Patterns::Bool(),
      "Enable simulations in a Galilean moving frame of reference");

    prm.enter_subsection("velocity");
    velocity->declare_parameters(prm, dim);
    prm.leave_subsection();

    prm.leave_subsection();
  }

  template <int dim>
  void
  ALE<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("ALE");
    enable = prm.get_bool("enable");

    prm.enter_subsection("velocity");
    velocity->parse_parameters(prm);
    prm.leave_subsection();

    prm.leave_subsection();
  }
} // namespace Parameters



#endif
