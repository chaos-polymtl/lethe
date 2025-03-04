// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef dealii_contact_force_h
#define dealii_contact_force_h

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

#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_contact_force.h>


// Tests (with common definitions)
#include <../tests/tests.h>

template <int dim, typename PropertiesIndex>
void 
set_particle_properties(int type,
                        Particles::ParticleIterator<dim> pit,
                        double particle_diameter,
                        Tensor<1, dim> v,
                        Tensor<1, dim> omega,
                        double mass)
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
















#endif //dealii_contact_force_h