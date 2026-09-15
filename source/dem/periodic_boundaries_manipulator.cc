// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/periodic_boundaries_manipulator.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <array>
#include <ranges>
#include <unordered_map>

using namespace dealii;

template <int dim>
PeriodicBoundariesManipulator<dim>::PeriodicBoundariesManipulator()
  : periodic_boundaries_enabled(false)
{}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::set_periodic_boundaries_information(
  const Parameters::PeriodicBoundaries &periodic_boundaries)
{
  // If this function is reached and the map is not empty
  if (periodic_boundaries.empty())
    return;

  // compute_periodic_offset_per_direction() sums every periodic boundary
  // pair's offset directly into a single Tensor<1, dim>, one component per
  // direction. This requires at most one periodic boundary pair per
  // direction: two pairs sharing a direction would silently add into the
  // same component and corrupt the periodic contact search.
  std::array<bool, dim> direction_is_used{};
  for (const auto &[id, boundary] : periodic_boundaries)
    {
      AssertThrow(boundary.direction < static_cast<unsigned int>(dim),
                  ExcMessage(
                    "PeriodicBoundariesManipulator: periodic direction " +
                    std::to_string(boundary.direction) +
                    " configured for boundary id " + std::to_string(id) +
                    " is out of range for a " + std::to_string(dim) +
                    "D simulation."));
      AssertThrow(
        !direction_is_used[boundary.direction],
        ExcMessage(
          "PeriodicBoundariesManipulator: more than one periodic boundary "
          "pair is configured for direction " +
          std::to_string(boundary.direction) +
          ". At most one periodic boundary pair per direction is "
          "supported."));
      direction_is_used[boundary.direction] = true;
    }

  periodic_boundaries_enabled = true;

  // Communicate to the action manager that there are periodic boundaries
  DEMActionManager::get_action_manager()->set_periodic_boundaries_enabled();

  this->periodic_boundaries = periodic_boundaries;

  // Initialize offset map
  this->periodic_offsets.clear();
}

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
  for (const auto &primary_mesh_id : periodic_boundaries | std::views::keys)
    offset_calculated[primary_mesh_id] = false;

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
                  const types::boundary_id face_boundary_id =
                    face->boundary_id();

                  // Check if the face matches any of the principal PB IDs
                  auto it = periodic_boundaries.find(face_boundary_id);
                  if (it != periodic_boundaries.end())
                    {
                      // Get direction corresponding to this periodic boundary
                      const unsigned int current_direction =
                        it->second.direction;

                      periodic_boundaries_cells_info_struct<dim>
                                   boundaries_information;
                      unsigned int face_id = cell->face_iterator_to_index(face);

                      get_periodic_boundaries_info(cell,
                                                   face_id,
                                                   boundaries_information);

                      // Insert into multimap
                      // One entry per boundary found
                      periodic_boundaries_cells_information.insert(
                        {boundaries_information.cell
                           ->global_active_cell_index(),
                         boundaries_information});

                      // Calculate offset if it has not yet been done for this
                      // periodic boundary id
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

  // Once periodic offsets calculated, combine them into a single
  // per-direction tensor
  this->compute_periodic_offset_per_direction();
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
  // Retrieve the correct offset for this specific boundary interaction.
  // boundaries_cells_content contains a single cell on a periodic boundary
  // and its associated periodic cell. The offset between these cells is
  // obtained from periodic_offsets to displace particles across these cells.
  // relevant_offset is a vector pointing from pb0 to pb1, so we check on which
  // side the current particles are and offset their positions accordingly
  // (either +offset or -offset).

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
      // return;
      AssertThrow(false,
                  ExcMessage(
                    "PeriodicBoundariesManipulator: no periodic offset found "
                    "for boundary id " +
                    std::to_string(boundaries_cells_content.boundary_id) +
                    ". This means map_periodic_cells() did not populate "
                    "periodic_offsets for this boundary."));
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
void
PeriodicBoundariesManipulator<dim>::compute_periodic_offset_per_direction()
{
  Tensor<1, dim> offset_per_direction;

  // Each entry of periodic_offsets is nonzero in exactly one component (the
  // direction of that periodic boundary pair, see get_periodic_boundaries_info
  // and map_periodic_cells), and set_periodic_boundaries_information()
  // guarantees at most one pair per direction, so summing them gives, for
  // every direction, the signed period of the domain along that direction
  // (0 if the direction is not periodic).
  for (auto const &[id, offset] : this->periodic_offsets)
    offset_per_direction += offset;

  this->periodic_offset_per_direction = offset_per_direction;
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
      // Iterate over all entries. Multimap correctly handles periodic corner
      // cells in higher dimensions. In 2D, a periodic corner cell would be
      // iterated over twice:  once for x-boundary, once for y-boundary.
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

          // Check and execute displacement of particles crossing a periodic
          // wall.
          // Case when particle goes from periodic boundary 0 to 1
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

          // Case when particle goes from periodic boundary 1 to 0
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
