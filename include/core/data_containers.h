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

#include <dem/particle_particle_contact_info_struct.h>
#include <dem/particle_point_line_contact_info_struct.h>
#include <dem/particle_wall_contact_info_struct.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_iterator.h>

#include <boost/range/adaptor/map.hpp>

#include <iostream>
#include <vector>

#ifndef lethe_data_containers_h
#  define lethe_data_containers_h

/**
 * This namespace defines all the data containers and required operations on
 * these containers for the DEM solver. At the moment, only overloaded operators
 * are here, but eventually all the data containers of the DEM solver will be
 * moved here
 */
namespace dem_data_containers
{
  /**
   * Operator overloading to enable using triangulation cells as map keys.
   */
  using namespace dealii;

  template <int dim>
  struct cell_comparison
  {
    bool
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
    bool
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
    bool
    operator()(
      const typename Triangulation<dim - 1, dim>::active_cell_iterator &cell_1,
      const typename Triangulation<dim - 1, dim>::active_cell_iterator &cell_2)
      const
    {
      return cell_1->global_active_cell_index() <
             cell_2->global_active_cell_index();
    }
  };

  template <int dim>
  struct dem_data_structures
  {
    typedef std::unordered_map<types::particle_index,
                               Particles::ParticleIterator<dim>>
      particle_index_iterator_map;

    typedef std::unordered_map<types::particle_index,
                               particle_point_line_contact_info_struct<dim>>
      particle_point_line_contact_info;

    typedef std::unordered_map<
      types::particle_index,
      std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
      particle_line_candidates;

    typedef std::unordered_map<
      types::particle_index,
      std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
      particle_point_candidates;

    typedef std::unordered_map<
      types::particle_index,
      std::unordered_map<types::boundary_id, Particles::ParticleIterator<dim>>>
      particle_floating_wall_candidates;

    typedef std::unordered_map<
      types::particle_index,
      std::unordered_map<types::boundary_id,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>,
                                    types::boundary_id,
                                    types::global_cell_index>>>
      particle_wall_candidates;

    typedef std::unordered_map<
      types::particle_index,
      std::map<types::boundary_id, particle_wall_contact_info_struct<dim>>>
      particle_wall_in_contact;

    typedef std::unordered_map<
      types::particle_index,
      std::unordered_map<types::particle_index,
                         particle_particle_contact_info_struct<dim>>>
      adjacent_particle_pairs;

    typedef std::vector<
      std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               std::unordered_map<types::particle_index,
                                  Particles::ParticleIterator<dim>>,
               cut_cell_comparison<dim>>>
      particle_floating_mesh_candidates;

    typedef std::vector<
      std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               std::unordered_map<types::particle_index,
                                  particle_wall_contact_info_struct<dim>>,
               cut_cell_comparison<dim>>>
      particle_floating_mesh_in_contact;

    typedef std::unordered_map<types::particle_index,
                               std::vector<types::particle_index>>
      particle_particle_candidates;

    typedef std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_neighbor_list;

    typedef std::vector<std::vector<
      std::pair<typename Triangulation<dim>::active_cell_iterator,
                typename Triangulation<dim - 1, dim>::active_cell_iterator>>>
      floating_mesh_information;

    typedef std::unordered_map<
      types::global_cell_index,
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_total_neighbor_list;

    typedef std::map<types::boundary_id, std::pair<Tensor<1, 3>, Point<3>>>
      boundary_points_and_normal_vectors;

    typedef std::map<unsigned int, std::map<types::boundary_id, Tensor<1, 3>>>
      vector_on_boundary;
  };

} // namespace dem_data_containers
#endif
