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
class RPTNodalReconstruction
{
public:
  RPTNodalReconstruction(std::vector<Detector<dim>> &detectors,
                         RPTCalculatingParameters   &rpt_parameters);

  void
  execute_nodal_reconstruction();

  void
  create_grid();

  void
  find_unknown_position();

  std::vector<Point<dim>>
  get_positions(
    IteratorRange<TriaIterator<CellAccessor<dim, dim>>> &cell_iterators,
    std::vector<int> parent_cell_indexes = {-1});

  std::vector<int>
  find_closer_cell(
    IteratorRange<TriaIterator<CellAccessor<dim, dim>>> &cell_iterators,
    std::vector<int> parent_cell_indexes = {-1});

  void
  calculate_counts(std::map<unsigned int, std::vector<double>> &index_counts);


  // void
  // output_results();

private:
  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler;
  FESystem<dim>      fe;
  PVDHandler         pvd_handler;
  MappingQ<dim>      mapping;
  QGauss<dim>        cell_quadrature;
  QGauss<dim - 1>    face_quadrature;

  std::vector<Point<dim>> positions;
  std::vector<Point<dim>> positions_coarse_mesh;

  std::map<unsigned int, std::vector<double>> map_vertices_index_coarse_mesh;
  RPTCalculatingParameters                    parameters;
  Parameters::RPTReconstructionParameters     reconstruction_parameters;
  std::vector<RadioParticle<dim>>             particle_positions;
  std::vector<Detector<dim>>                  detectors;
  Vector<double>                              map_counts;
  std::vector<double>                         reconstruction_counts;
  std::map<unsigned int, std::vector<double>> map_vertices_index;
};



#endif // LETHE_MAP_H
