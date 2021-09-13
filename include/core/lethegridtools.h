/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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
 * -------------------------------------------------------------------*/

#ifndef lethe_lethegridtools_h
#define lethe_lethegridtools_h


#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

#include <unordered_set>

using namespace dealii;
namespace LetheGridTools
{
    /**
     * @brief
     * Map the vertex index to the cells that includes that vertex.
     *
     * @param dof_handler DofHandler of the triangulation on which to create the map
     *
     * @param vertices_cell_map Object to which the data will be written
     */

    template <int dim>
    void
    vertices_cell_mapping(const DoFHandler<dim> &dof_handler,
            std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> & vertices_cell_map);

    /**
     * @brief
     * Return the cell around a point by starting from the least refined level and iterating over children
     *
     * @param cell The initial cell. We suspect the point of being in one of the neighbours of this cell.
     *
     * @param point The point that we want to find the cell that contains it
     */
    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                     Point<dim>             &point);

    /**
     * @brief
     * Return the cell around a point based on a initial guess of a closed cell
     * (look in the neighbors of this cell)
     *
     * @param cell The initial cell. We suspect the point of being in one of the neighbours of this cell.
     *
     * @param point The point that we want to find the cell that contains it
     */

    template<int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree_with_guess(const DoFHandler<dim> &dof_handler,
                                                const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                const Point<dim> &point);

    template<int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_neighbors(const DoFHandler<dim> &dof_handler,
                                          std::map<unsigned int, std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          const Point<dim> &point);

    /**
      * @brief
      *Return a vector of cells around a cell including vertex neighbors
      *
      * @param cell The initial cell. we want to know all the cells that share a vertex with this cell.
      */
    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_cell(std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
            const typename DoFHandler<dim>::active_cell_iterator &cell);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_flat_cell(const DoFHandler<dim> &                                   dof_handler,
                                const typename DoFHandler<dim - 1,dim>::active_cell_iterator &cell,
                                std::map<unsigned int,
                                        std::set<typename DoFHandler<dim>::active_cell_iterator>>
                                &vertices_cell_map);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_edge(const DoFHandler<dim> &dof_handler,
                           std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                           const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                           Point<dim> &point_2);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_in_cells(const DoFHandler<dim> &dof_handler_1,
            const typename DoFHandler<dim>::active_cell_iterator &cell);

    template <int dim>
    bool
    cell_cut_by_flat(const typename DoFHandler<dim>::active_cell_iterator &cell,const typename DoFHandler<dim-1,dim>::active_cell_iterator &cell_flat);

    template <int dim>
    bool
    cell_pierced_by_edge(const typename DoFHandler<dim>::active_cell_iterator &cell, const TriaIterator<CellAccessor<1, dim>> &cell_edge);

    template <int dim>
    bool
    cell_pierced_by_edge(const typename DoFHandler<dim>::active_cell_iterator &cell, Point<dim> point_1, Point<dim> point_2);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    move_grid(Triangulation<dim> mesh,Tensor<2,dim+1> displacement);

    template  <int spacedim, int structdim>
    std::pair<std::pair<Point<spacedim>,bool>,Tensor<1,spacedim>>
    project_to_d_linear_object(const typename DoFHandler<structdim,spacedim>::active_cell_iterator &       object,
                               const Point<spacedim> &trial_point);

    template  <int spacedim, int structdim>
    double
    dist_based_on_reference_point(const typename DoFHandler<structdim,spacedim>::active_cell_iterator&      object,
                                        const Point<spacedim> &trial_point, const Point<structdim> &xi);

    template < int spacedim, int structdim>
    Tensor<1,structdim>
    grad_dist_based_on_reference_point(const typename DoFHandler<structdim,spacedim>::active_cell_iterator&      object,
                                                       const Point<spacedim> &trial_point, const Point<structdim> &xi);


    template <int dim>
    struct hash_cell
    {
      std::size_t
      operator()(const typename DoFHandler<dim>::active_cell_iterator &cell)
        const noexcept
      {
        return  std::hash<int>()(cell->global_active_cell_index());
      }
    };

    template <int dim>
    struct equal_cell
    {
      std::size_t
      operator()(const typename DoFHandler<dim>::active_cell_iterator &cell1,
                 const typename DoFHandler<dim>::active_cell_iterator &cell2)
        const noexcept
      {
        return cell1->global_active_cell_index()==cell2->global_active_cell_index();
      }
    };

}


#endif //lethe_lethegridtools_h

