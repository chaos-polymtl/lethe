/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 */

#include <dem/boundary_cells_info_struct.h>
#include <dem/particle_particle_contact_info.h>
#include <dem/particle_point_line_contact_info_struct.h>
#include <dem/particle_wall_contact_info.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_iterator.h>

#include <boost/range/adaptor/map.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
#include <ankerl/unordered_dense.h>
#pragma GCC diagnostic pop

#include <iostream>
#include <unordered_map>
#include <vector>

#ifndef lethe_data_containers_h
#  define lethe_data_containers_h

/**
 * This namespace defines all the data containers and required operations on
 * these containers for the DEM solver. At the moment, only overloaded operators
 * are here, but eventually all the data containers of the DEM solver will be
 * moved here
 */
namespace DEM
{
  /**
   * Defined global face id as a type to ensure that we know what the
   * structures use as index. This is more verbose than unsigned int
   */
  typedef unsigned int global_face_id;


  /**
   * Operator overloading to enable using triangulation cells as map keys.
   */
  using namespace dealii;

  template <int dim>
  struct cell_comparison
  {
    inline bool
    operator()(
      const typename Triangulation<dim>::active_cell_iterator &cell_1,
      const typename Triangulation<dim>::active_cell_iterator &cell_2) const
    {
      return cell_1->global_active_cell_index() <
             cell_2->global_active_cell_index();
    }
  };

  /**
   * Operator overloading to enable using particle iterators as map keys.
   */

  template <int dim>
  struct particle_comparison
  {
    inline bool
    operator()(const Particles::ParticleIterator<dim> &particle_1,
               const Particles::ParticleIterator<dim> &particle_2)
    {
      return particle_1->id() < particle_2->id();
    }
  };

  /**
   * Operator overloading to enable using triangulation cells as map keys.
   */
  using namespace dealii;

  template <int dim>
  struct cut_cell_comparison
  {
    inline bool
    operator()(
      const typename Triangulation<dim - 1, dim>::active_cell_iterator &cell_1,
      const typename Triangulation<dim - 1, dim>::active_cell_iterator &cell_2)
      const
    {
      return cell_1->global_active_cell_index() <
             cell_2->global_active_cell_index();
    }
  };

  // DEM data structure containers
  // To improve readability of containers, there's a simple description above
  // all declarations which follow these representation :
  // <> : map or unordered map
  // () : tuple or pair
  // [] : vector
  template <int dim>
  struct dem_data_structures
  {
    // <particle id, particle iterator>
    typedef ankerl::unordered_dense::map<types::particle_index,
                                         Particles::ParticleIterator<dim>>
      particle_index_iterator_map;

    // <particle id, point-line info>
    typedef std::unordered_map<types::particle_index,
                               particle_point_line_contact_info_struct<dim>>
      particle_point_line_contact_info;

    // <particle id, (particle iterator, point, point)>
    typedef std::unordered_map<
      types::particle_index,
      std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
      particle_line_candidates;

    // <particle id, (particle iterator, point)>
    typedef std::unordered_map<
      types::particle_index,
      std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
      particle_point_candidates;

    // <particle id, <global face id, particle iterator>>
    typedef ankerl::unordered_dense::map<
      types::particle_index,
      std::unordered_map<global_face_id, Particles::ParticleIterator<dim>>>
      particle_floating_wall_candidates;

    // <particle id, <global_face_id, (particle iterator, tensor, point,
    // boundary id, cell id)>>
    typedef ankerl::unordered_dense::map<
      types::particle_index,
      std::unordered_map<global_face_id,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>,
                                    global_face_id>>>
      particle_wall_candidates;

    // <particle id, <global_face_id, particle-wall info>>
    typedef ankerl::unordered_dense::map<
      types::particle_index,
      std::unordered_map<global_face_id, particle_wall_contact_info<dim>>>
      particle_wall_in_contact;

    // <particle id, <particle id, particle-particle info>>
    typedef ankerl::unordered_dense::map<
      types::particle_index,
      ankerl::unordered_dense::map<types::particle_index,
                                   particle_particle_contact_info<dim>>>
      adjacent_particle_pairs;

    // <cell iterator, <particle id, particle iterator>>
    typedef std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
                     std::unordered_map<types::particle_index,
                                        Particles::ParticleIterator<dim>>,
                     /* mapped_type */ cut_cell_comparison<dim>>
      particle_floating_wall_from_mesh_candidates;

    // <cell iterator, <particle id, particle-wall info>>
    typedef std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
                     std::unordered_map<types::particle_index,
                                        particle_wall_contact_info<dim>>,
                     /* mapped_type */ cut_cell_comparison<dim>>
      particle_floating_wall_from_mesh_in_contact;

    // [<cell iterator, <particle id, particle iterator>>]
    typedef std::vector<particle_floating_wall_from_mesh_candidates>
      particle_floating_mesh_candidates;

    // [<cell iterator, <particle id, particle-wall info>>]
    typedef std::vector<particle_floating_wall_from_mesh_in_contact>
      particle_floating_mesh_in_contact;

    // <particle id, [particle id]>
    typedef ankerl::unordered_dense::map<types::particle_index,
                                         std::vector<types::particle_index>>
      particle_particle_candidates;

    // [[cell iterator]]
    typedef std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_neighbor_list;

    // [[(cell iterator, cell iterator)]]
    typedef std::vector<std::vector<
      std::pair<typename Triangulation<dim>::active_cell_iterator,
                typename Triangulation<dim - 1, dim>::active_cell_iterator>>>
      floating_mesh_information;

    // <cell id, [cell iterator]>
    typedef ankerl::unordered_dense::map<
      types::global_cell_index,
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_total_neighbor_list;

    // <global_face_id, (tensor, point)>
    typedef ankerl::unordered_dense::map<global_face_id,
                                         std::pair<Tensor<1, 3>, Point<3>>>
      boundary_points_and_normal_vectors;

    // <unsigned int, <global_face_id, tensor>>
    typedef ankerl::unordered_dense::
      map<unsigned int, std::map<types::boundary_id, Tensor<1, 3>>>
        vector_on_boundary;

    // [cell iterators]
    typedef std::vector<typename Triangulation<dim>::active_cell_iterator>
      cell_vector;

    // [cell iterators]
    typedef std::set<typename Triangulation<dim>::active_cell_iterator>
      cell_set;

    // <cell id, periodic cells info>
    typedef ankerl::unordered_dense::
      map<types::global_cell_index, periodic_boundaries_cells_info_struct<dim>>
        periodic_boundaries_cells_info;

    // <cell id, integer value>
    typedef ankerl::unordered_dense::map<types::global_cell_index, unsigned int>
      cell_index_int_map;
  };

} // namespace DEM
#endif
