// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_contact_info_h
#define lethe_contact_info_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

using namespace dealii;

/**
 * @brief Handle the information related to the calculation of the
 * particle-particle contact force. Notably it is responsible for storing
 * information that has to be preserved over multiple iterations of a contact,
 * namely everything related to tangential displacements
 */
template <int dim>
struct particle_particle_contact_info
{
  Particles::ParticleIterator<dim> particle_one;
  Particles::ParticleIterator<dim> particle_two;
  Tensor<1, 3>                     tangential_displacement;
  Tensor<1, 3>                     rolling_resistance_spring_torque;
};

/**
 * @brief Handle the information related to the calculation of the
 * particle-wall contact forces.
 */
template <int dim>
struct particle_wall_contact_info
{
  Particles::ParticleIterator<dim> particle;
  Tensor<1, 3>                     normal_vector;
  Point<3>                         point_on_boundary;
  types::boundary_id               boundary_id;
  Tensor<1, 3>                     tangential_displacement;
  Tensor<1, 3>                     rolling_resistance_spring_torque;
};

/**
 * @brief Handle information related to the calculation of the particle-line
 * contact forces.
 */
template <int dim>
struct particle_line_contact_info
{
  Particles::ParticleIterator<dim> particle;
  Point<3>                         point_one;
  Point<3>                         point_two;
};

/**
 * @brief Handle information related to the calculation of the particle-point
 * contact forces.
 */
template <int dim>
struct particle_point_contact_info
{
  Particles::ParticleIterator<dim> particle;
  Point<3>                         point;
};

/**
 * @brief Handle information related to the cell-line matching.
 */
template <int dim>
struct cell_line_info
{
  typename Triangulation<dim>::active_cell_iterator cell;
  Point<3>                                          point_one;
  Point<3>                                          point_two;
};

/**
 * @brief Handle information related to the cell-point matching.
 */
template <int dim>
struct cell_point_info
{
  typename Triangulation<dim>::active_cell_iterator cell;
  Point<3>                                          point;
};

#endif
