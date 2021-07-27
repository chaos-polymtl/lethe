/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Authors: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

#ifndef lethe_rpt_map_h
#define lethe_rpt_map_h

#include <deal.II/base/config.h>

#include <core/pvd_handler.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/tria.h>

#include <rpt/detector.h>
#include <rpt/parameters_rpt.h>
#include <rpt/radioactive_particle.h>
#include <rpt/rpt_calculating_parameters.h>

#include <vector>

using namespace dealii;

template <int dim>
class RPTMap
{
public:
  RPTMap(std::vector<Detector<dim>> &detectors,
         RPTCalculatingParameters &  rpt_parameters);

  void
  add_calculated_counts(std::vector<double> calculated_counts);

  void
  create_grid();

  std::vector<Point<dim>>
  get_positions();

  void
  output_results();

private:
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler;
  FESystem<dim>      fe;
  PVDHandler         pvd_handler;
  MappingQ<dim>      mapping;
  QGauss<dim>        cell_quadrature;
  QGauss<dim - 1>    face_quadrature;

  std::vector<Point<dim>>                 positions;
  Parameters::RPTParameters               parameters;
  Parameters::RPTReconstructionParameters reconstruction_parameters;
  std::vector<RadioParticle<dim>>         particle_positions;
  std::vector<Detector<dim>>              detectors_positions;
  Vector<double>                          counts;
  unsigned int                            p = 1;
};



#endif // LETHE_MAP_H
