// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_data_containers_h
#define lethe_data_containers_h

#include <core/utilities.h>

#include <dem/boundary_cells_info_struct.h>
#include <dem/contact_info.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_iterator.h>

#include <unordered_map>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
#include <ankerl/unordered_dense.h>
#pragma GCC diagnostic pop

/**
 * @brief This namespace defines all the data containers and required operations on
 * these containers for the DEM solver. At the moment, only overloaded operators
 * are here, but eventually all the data containers of the DEM solver will be
 * moved here
 */
namespace DEM
{
  /**
   * @brief Defined global face id as a type to ensure that we know what the
   * structures use as index. This is more verbose than unsigned int
   */
  using global_face_id = unsigned int;

  using namespace dealii;



  /**
   * @brief DEM data structure containers
   * To improve readability of containers, there's a simple description above
   * all declarations which follow these representation :
   * <> : map or unordered map or unordered multimap
   * () : tuple or pair
   * [] : vector
   */
  template <int dim>
  struct dem_data_structures
  {
    // <particle id, particle iterator>
    using particle_index_iterator_map =
      ankerl::unordered_dense::map<types::particle_index,
                                   Particles::ParticleIterator<dim>>;

    // <particle id, line info>
    using particle_line_in_contact =
      std::unordered_map<types::particle_index,
                         particle_line_contact_info<dim>>;

    // <particle id, point info>
    using particle_point_in_contact =
      std::unordered_map<types::particle_index,
                         particle_point_contact_info<dim>>;

    // <particle id, line info>
    using particle_line_candidates =
      std::unordered_map<types::particle_index,
                         particle_line_contact_info<dim>>;

    // <particle id, (particle iterator, point)>
    using particle_point_candidates =
      std::unordered_map<types::particle_index,
                         particle_point_contact_info<dim>>;

    // <particle id, <global face id, particle iterator>>
    using particle_floating_wall_candidates = ankerl::unordered_dense::map<
      types::particle_index,
      std::unordered_map<global_face_id, Particles::ParticleIterator<dim>>>;

    // <particle id, <global_face_id, (particle iterator, tensor, point,
    // boundary id, cell id)>>
    using particle_wall_candidates = ankerl::unordered_dense::map<
      types::particle_index,
      ankerl::unordered_dense::map<global_face_id,
                                   std::tuple<Particles::ParticleIterator<dim>,
                                              Tensor<1, dim>,
                                              Point<dim>,
                                              global_face_id>>>;

    // <particle id, <global_face_id, particle-wall info>>
    using particle_wall_in_contact = ankerl::unordered_dense::map<
      types::particle_index,
      ankerl::unordered_dense::map<global_face_id,
                                   particle_wall_contact_info<dim>>>;

    // <particle id, particle-particle info>
    using particle_contact_info =
      ankerl::unordered_dense::map<types::particle_index,
                                   particle_particle_contact_info<dim>>;

    // <particle id, periodic particle-particle info>
    using periodic_particle_contact_info = ankerl::unordered_dense::
      map<types::particle_index, periodic_particle_particle_contact_info<dim>>;

    // Value stored per "particle one" in the adjacency containers below. The
    // iterator to particle one is shared by all its contacts, so it is stored
    // once here instead of being duplicated in every particle-particle contact
    // info. second_particles maps uses each particle two id (key)
    // to access the contact info (value) between particle one and two.
    // (particle one iterator, <particle id, particle-particle info>)
    struct adjacent_particle_pairs_value
    {
      Particles::ParticleIterator<dim> particle_one;
      particle_contact_info            second_particles;
    };

    // Periodic counterpart of adjacent_particle_pairs_value.
    // (particle one iterator, <particle id, periodic particle-particle info>)
    struct periodic_adjacent_particle_pairs_value
    {
      Particles::ParticleIterator<dim> particle_one;
      periodic_particle_contact_info   second_particles;
    };

    // <particle id, (particle one iterator, <particle id, particle-particle
    // info>)>
    using adjacent_particle_pairs =
      ankerl::unordered_dense::map<types::particle_index,
                                   adjacent_particle_pairs_value>;

    // <particle id, (particle one iterator, <particle id, periodic
    // particle-particle info>)>
    using periodic_adjacent_particle_pairs =
      ankerl::unordered_dense::map<types::particle_index,
                                   periodic_adjacent_particle_pairs_value>;

    // <cell iterator, <particle id, particle iterator>>
    using particle_floating_wall_from_mesh_candidates =
      std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               ankerl::unordered_dense::map<types::particle_index,
                                            Particles::ParticleIterator<dim>>,
               cut_cell_comparison<dim>>;

    // <cell iterator, <particle id, particle-wall info>>
    using particle_triangle_cell_from_mesh_potentially_in_contact =
      std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               ankerl::unordered_dense::map<types::particle_index,
                                            particle_wall_contact_info<dim>>,
               cut_cell_comparison<dim>>;

    // [<cell iterator, <particle id, particle iterator>>]
    using particle_floating_mesh_candidates =
      std::vector<particle_floating_wall_from_mesh_candidates>;

    // [cell iterator]
    using triangulation_cell_vector =
      std::vector<typename Triangulation<dim - 1, dim>::active_cell_iterator>;

    // <cell iterator, [cell iterator]>
    using triangulation_cell_map =
      std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               triangulation_cell_vector,
               cut_cell_comparison<dim>>;

    // [<cell iterator, <particle id, particle-wall info>>]
    using particle_floating_mesh_potentially_in_contact =
      std::vector<particle_triangle_cell_from_mesh_potentially_in_contact>;

    // <particle id, [particle id]>
    using particle_particle_candidates =
      ankerl::unordered_dense::map<types::particle_index,
                                   std::vector<types::particle_index>>;

    // [[cell iterator]]
    using cells_neighbor_list = std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>>;

    // [[(background cell iterator, solid surface cell iterator)]]
    using solid_surfaces_mesh_information = std::vector<std::vector<
      std::pair<typename Triangulation<dim>::active_cell_iterator,
                typename Triangulation<dim - 1, dim>::active_cell_iterator>>>;

    // [[(background cell iterator, solid volume cell iterator)]]
    using solid_volumes_mesh_info = std::vector<std::vector<
      std::pair<typename Triangulation<dim>::active_cell_iterator,
                typename Triangulation<dim, dim>::active_cell_iterator>>>;

    // <cell id, [cell iterator]>
    using cells_total_neighbor_list = ankerl::unordered_dense::map<
      types::global_cell_index,
      std::vector<typename Triangulation<dim>::active_cell_iterator>>;

    // <global_face_id, (tensor, point)>
    using boundary_points_and_normal_vectors =
      ankerl::unordered_dense::map<global_face_id,
                                   std::pair<Tensor<1, 3>, Point<3>>>;

    // <unsigned int, <global_face_id, tensor>>
    using vector_on_boundary =
      ankerl::unordered_dense::map<unsigned int,
                                   std::map<types::boundary_id, Tensor<1, 3>>>;

    // [cell iterators]
    using cell_vector =
      std::vector<typename Triangulation<dim>::active_cell_iterator>;

    // [cell iterators]
    using cell_set =
      std::set<typename Triangulation<dim>::active_cell_iterator>;

    // <cell id, periodic cells info>
    using periodic_boundaries_cells_info =
      std::unordered_multimap<types::global_cell_index,
                              periodic_boundaries_cells_info_struct<dim>>;

    // <cell id, integer value>
    using cell_index_int_map =
      ankerl::unordered_dense::map<types::global_cell_index, unsigned int>;
  };

} // namespace DEM
#endif
