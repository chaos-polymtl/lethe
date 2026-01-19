// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef particle_handler_conversion_h
#define particle_handler_conversion_h


#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

using namespace dealii;

/** @brief Converts an input particle handler with a given set of input properties
 * into an output particle handler with a different set of output properties.
 * This is used, notably, to convert DEM simulation particle_handlers into
 * CFD-DEM particle_handlers.
 *
 * @tparam dim the number of dimensions
 * @tparam input the properties index used for the ph_in
 * @tparam output the properties index used for the ph_out
 *
 * @param[in] triangulation The triangulation on which the ph_out is stored
 * @param[in] ph_in The input particle_handler that has properties identified by
 * the input.
 * @param[out] ph_out The output particle_handler that has properties identified
 *by the output.
 */
template <int dim, typename input, typename output>
void
convert_particle_handler(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim>           &ph_in,
  Particles::ParticleHandler<dim>                 &ph_out)
{
  // Pre-allocate the vector of particle properties and location of ph_out
  std::vector<std::vector<double>> ph_out_properties;
  std::vector<Point<dim>>          ph_out_points;
  ph_out_properties.reserve(ph_in.n_locally_owned_particles());
  ph_out_points.reserve(ph_in.n_locally_owned_particles());

  // We loop through all particles of the input particle handler and gather the
  // properties and locations
  for (auto particle = ph_in.begin(); particle != ph_in.end(); ++particle)
    {
      // Get and store particle location
      auto particle_location = particle->get_location();
      ph_out_points.emplace_back(particle_location);

      // Get the particle properties
      auto particle_properties = particle->get_properties();
      // Give a default value of zero to all properties. This will account for
      // increasing number of properties
      std::vector<double> ph_out_particle_properties(output::n_properties, 0.);
      ph_out_particle_properties[output::type] =
        particle_properties[input::type];
      ph_out_particle_properties[output::dp] = particle_properties[input::dp];

      ph_out_particle_properties[output::v_x] = particle_properties[input::v_x];
      ph_out_particle_properties[output::v_y] = particle_properties[input::v_y];
      ph_out_particle_properties[output::v_z] = particle_properties[input::v_z];

      ph_out_particle_properties[output::omega_x] =
        particle_properties[input::omega_x];
      ph_out_particle_properties[output::omega_y] =
        particle_properties[input::omega_y];
      ph_out_particle_properties[output::omega_z] =
        particle_properties[input::omega_z];

      ph_out_particle_properties[output::mass] =
        particle_properties[input::mass];

      ph_out_properties.emplace_back(ph_out_particle_properties);
    }


  // We use the insert global functionalities to insert the particles in the
  // output particle handler. This is not optimal since it will recalculate the
  // reference location of every particle, but since this operation is only done
  // once per simulation, this cost should be negligible.

  // Calculate global bounding box for global insertion
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(triangulation.get_mpi_communicator(),
                               my_bounding_box);

  // Do the global insertion
  ph_out.insert_global_particles(ph_out_points,
                                 global_bounding_boxes,
                                 ph_out_properties);
}


#endif
