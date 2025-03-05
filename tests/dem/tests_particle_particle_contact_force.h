// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef tests_particle_particle_contact_force_h
#define tests_particle_particle_contact_force_h


// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/particle_particle_contact_force.h>


// Tests (with common definitions)
#include <../tests/tests.h>


/**
 * @brief Return the particle iterator with the inserted particle, according to the position of the particle.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param triangulation Triangulation to access the information of the cells.
 * @param position Postion of the particle.
 * @param id Id of the particle.
 */
template <int dim>
Particles::ParticleIterator<dim>
construct_particle_iterator(Particles::ParticleHandler<dim> &particle_handler,
                            parallel::distributed::Triangulation<dim> &triangulation,
                            Point<3> &position,
                            int id)
{
  typename Triangulation<dim>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(triangulation, position);
  Particles::Particle<dim> particle(position, position, id);
  return particle_handler.insert_particle(particle, cell);
}


/**
 * @brief Set the properties of the particle in the PropertiesIndex.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param pit Particle Iterator where the particle is inserted.
 * @param type Particle type.
 * @param particle_diameter Particle diameter.
 * @param mass Particle mass.
 * @param v Particle velocity.
 * @param omega Particle angular velocity.
 */
template <int dim, typename PropertiesIndex>
void 
set_particle_properties(Particles::ParticleIterator<dim> &pit,
                        int type,
                        double particle_diameter,
                        double mass,
                        Tensor<1, dim> &v,
                        Tensor<1, dim> &omega)
{
  pit->get_properties()[PropertiesIndex::type]    = type;
  pit->get_properties()[PropertiesIndex::dp]      = particle_diameter;
  pit->get_properties()[PropertiesIndex::v_x]     = v[0];
  pit->get_properties()[PropertiesIndex::omega_x] = omega[0];
  if (dim>1) {pit->get_properties()[PropertiesIndex::v_y]     = v[1];
              pit->get_properties()[PropertiesIndex::omega_y] = omega[1];}
  if (dim>2) {pit->get_properties()[PropertiesIndex::v_z]     = v[2];
              pit->get_properties()[PropertiesIndex::omega_z] = omega[2];}
  pit->get_properties()[PropertiesIndex::mass]    = mass;

}


#endif //tests_particle_particle_contact_force_h