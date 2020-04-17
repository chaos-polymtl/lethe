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

#ifndef LETHE_MANIFOLDS_H
#define LETHE_MANIFOLDS_H

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria_base.h>


#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

using namespace dealii;

namespace Parameters
{
  class Manifolds
  {
  public:
    enum class ManifoldType
    {
      none,
      spherical,
      iges
    };

    // ID of boundary condition
    std::vector<unsigned int> id;

    // List of boundary type for each number
    std::vector<ManifoldType> types;

    // Arguments of manifold
    std::vector<double> arg1;
    std::vector<double> arg2;
    std::vector<double> arg3;
    std::vector<double> arg4;
    std::vector<double> arg5;
    std::vector<double> arg6;

    // File names for cad manifolds
    std::vector<std::string> cad_files;

    // Number of boundary conditions
    unsigned int size;
    unsigned int max_size;

    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);

    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace Parameters

void attach_cad_to_manifold(std::shared_ptr<parallel::DistributedTriangulationBase<2>> triangulation,std::string cad_name, int manifold_id);

void attach_cad_to_manifold(std::shared_ptr<parallel::DistributedTriangulationBase<3>> triangulation, std::string cad_name, int manifold_id);


#endif
