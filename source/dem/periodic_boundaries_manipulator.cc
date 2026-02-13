// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/periodic_boundaries_manipulator.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <unordered_map>

// TODO: data structure periodic_boundaries_cells_info must be defined
// as std::unordered_multimap

using namespace dealii;

template <int dim>
PeriodicBoundariesManipulator<dim>::PeriodicBoundariesManipulator()
  : periodic_boundaries_enabled(false)
{}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::get_periodic_boundaries_info(
  typename Triangulation<dim>::cell_iterator  cell,
  const unsigned int                          face_id,
  periodic_boundaries_cells_info_struct<dim> &boundaries_information)
{
  // Initialize a simple quadrature for the system. This will be used to
  // obtain a single sample point on the boundaries faces.
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);
  Tensor<1, dim>    normal_vector;
  Point<dim>        quad_point;

  // Store information of the cell
  boundaries_information.cell        = cell;
  boundaries_information.boundary_id = cell->face(face_id)->boundary_id();

  // Store information of the periodic cell
  boundaries_information.periodic_cell = cell->periodic_neighbor(face_id);

  // Find the normal vector of the boundary face and point
  fe_face_values.reinit(cell, face_id);
  boundaries_information.normal_vector = fe_face_values.normal_vector(0);
  boundaries_information.point_on_face = fe_face_values.quadrature_point(0);

  // Find the normal vector of the corresponding periodic boundary face and
  // point with the reinit FE face values
  fe_face_values.reinit(boundaries_information.periodic_cell,
                        cell->periodic_neighbor_face_no(face_id));
  normal_vector = fe_face_values.normal_vector(0);
  quad_point    = fe_face_values.quadrature_point(0);
  boundaries_information.periodic_normal_vector = normal_vector;
  boundaries_information.point_on_periodic_face = quad_point;
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::map_periodic_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information)
{
  if (!periodic_boundaries_enabled)
    return;

  periodic_boundaries_cells_information.clear();
  periodic_offsets.clear();

  // Temp storage to calculate offsets only once per ID
  std::map<types::boundary_id, bool> offset_calculated;
  for (auto [bc_index, id] : periodic_boundaries_ids)
    offset_calculated[id] = false;

  // Iterating over the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          if (cell->at_boundary())
            {
              // Iterating over cell faces
              for (const auto &face : cell->face_iterators())
                {
                  unsigned int face_boundary_id = face->boundary_id();

                  // Check if face matches any of the PB IDs
                  for (auto const &[bc_index, primary_mesh_id] :
                       periodic_boundaries_ids)
                    {
                      if (face_boundary_id == primary_mesh_id)
                        {
                          // Get direction corresponding to this BC index
                          unsigned int current_direction =
                            directions.at(bc_index);

                          periodic_boundaries_cells_info_struct<dim>
                            boundaries_information;
                          unsigned int face_id =
                            cell->face_iterator_to_index(face);

                          get_periodic_boundaries_info(cell,
                                                       face_id,
                                                       boundaries_information);

                          // Insert into multimap
                          // One entry per boundary found
                          periodic_boundaries_cells_information.insert(
                            {boundaries_information.cell
                               ->global_active_cell_index(),
                             boundaries_information});

                          // Calculate offset if not yet done for this PB ID
                          if (!offset_calculated[face_boundary_id])
                            {
                              Tensor<1, dim> offset;
                              offset[current_direction] =
                                boundaries_information
                                  .point_on_periodic_face[current_direction] -
                                boundaries_information
                                  .point_on_face[current_direction];

                              periodic_offsets[face_boundary_id]  = offset;
                              offset_calculated[face_boundary_id] = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::check_and_move_particles(
  const periodic_boundaries_cells_info_struct<dim> &boundaries_cells_content,
  const bool                                       &particles_in_pb0_cell,
  typename Particles::ParticleHandler<dim>::particle_iterator_range
       &particles_in_cell,
  bool &particle_has_been_moved)
{
  // Retrieve correct offset for this specific boundary interaction.
  // If particles are in pb0 cell, we use the stored offset map.
  // If they are in the periodic neighbor (pb1), we infer the ID from
  // the stored pb0 ID. The struct stores the 'boundary_id' of the face.

  Tensor<1, dim> relevant_offset;

  // boundaries_cells_content.boundary_id corresponds to the face on the
  // "main" side (pb0)
  auto offset_it = periodic_offsets.find(boundaries_cells_content.boundary_id);
  if (offset_it != periodic_offsets.end())
    {
      relevant_offset = offset_it->second;
    }
  else
    {
      // Fallback/error, map_periodic_cells should populate this
      return;
    }

  for (auto particle = particles_in_cell.begin();
       particle != particles_in_cell.end();
       ++particle)
    {
      // Get the current particle location
      Point<dim> particle_position = particle->get_location();

      // Initialize points and normal vector related of the cell that contains
      // the particle
      Point<dim>     point_on_face;
      Tensor<1, dim> normal_vector, distance_between_faces;

      if (particles_in_pb0_cell)
        {
          point_on_face          = boundaries_cells_content.point_on_face;
          normal_vector          = boundaries_cells_content.normal_vector;
          distance_between_faces = relevant_offset;
        }
      else
        {
          point_on_face = boundaries_cells_content.point_on_periodic_face;
          normal_vector = boundaries_cells_content.periodic_normal_vector;
          distance_between_faces = -relevant_offset;
        }

      // Calculate distance between particle position and the face of
      // boundary cell d = n•(pt_particle - pt_face)
      double distance_with_face =
        scalar_product(particle_position - point_on_face, normal_vector);

      // If distance >= 0, particle is outside of cell (or on face).
      // If so, particle location is modified to get moved into the periodic
      // cell.
      if (distance_with_face >= 0.0)
        {
          // Move particle outside the current cell to the periodic cell.
          particle_position += distance_between_faces;
          particle->set_location(particle_position);

          // Update flag to indicate that particle has been moved
          particle_has_been_moved = true;
        }
    }
}

template <int dim>
bool
PeriodicBoundariesManipulator<dim>::execute_particles_displacement(
  const Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information)
{
  if (!periodic_boundaries_enabled)
    return false;

  bool particle_has_been_moved = false;

  if (!periodic_boundaries_cells_information.empty())
    {
      // Iterate over all entries. If using multimap, this correctly handles
      // corner cells twice :once for x-boundary, once for y-boundary (in 2D)
      for (auto boundaries_cells_information_iterator =
             periodic_boundaries_cells_information.begin();
           boundaries_cells_information_iterator !=
           periodic_boundaries_cells_information.end();
           ++boundaries_cells_information_iterator)
        {
          // Get the cell and periodic cell from map
          auto boundaries_cells_content =
            boundaries_cells_information_iterator->second;
          auto cell          = boundaries_cells_content.cell;
          auto periodic_cell = boundaries_cells_content.periodic_cell;

          // Main cell (pb0 side)
          if (cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler.particles_in_cell(cell);
              const bool particles_exist_in_main_cell =
                !particles_in_cell.empty();

              if (particles_exist_in_main_cell)
                {
                  bool particles_in_pb0_cell = true;
                  check_and_move_particles(boundaries_cells_content,
                                           particles_in_pb0_cell,
                                           particles_in_cell,
                                           particle_has_been_moved);
                }
            }

          // Periodic neighbor (pb1 side)
          if (periodic_cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_cell =
                  particle_handler.particles_in_cell(periodic_cell);
              const bool particles_exist_in_periodic_cell =
                !particles_in_periodic_cell.empty();

              if (particles_exist_in_periodic_cell)
                {
                  bool particles_in_pb0_cell = false;
                  check_and_move_particles(boundaries_cells_content,
                                           particles_in_pb0_cell,
                                           particles_in_periodic_cell,
                                           particle_has_been_moved);
                }
            }
        }
    }

  return particle_has_been_moved;
}

template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;
