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
 */


#ifndef LETHE_PARTICLE_VISUALIZATION_H
#define LETHE_PARTICLE_VISUALIZATION_H

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include <rpt/rpt_calculating_parameters.h>

using namespace dealii;

template <int dim>
class ParticleVisualization
{
public:
  ParticleVisualization(Triangulation<dim> &     background_triangulation,
                        std::string &            filename,
                        std::vector<Point<dim>> &positions,
                        std::vector<std::vector<double>> &counts);

  void
  build_visualization_files();

  void
  output_particles(unsigned int it);



private:
  std::vector<Point<dim>>          particle_positions;
  std::vector<std::vector<double>> particle_counts;
  double                           dt;
  std::string                      visualization_filename;
  Particles::ParticleHandler<dim>  particle_handler;
  DoFHandler<dim>                  empty_dof_handler;
  MappingQGeneric<dim>             mapping;
};



#endif // LETHE_PARTICLE_VISUALIZATION_H
