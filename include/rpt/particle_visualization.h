// SPDX-FileCopyrightText: Copyright (c) 2021-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef LETHE_PARTICLE_VISUALIZATION_H
#define LETHE_PARTICLE_VISUALIZATION_H

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q_generic.h>

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
  ParticleVisualization(Triangulation<dim>      &background_triangulation,
                        std::string             &filename,
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
