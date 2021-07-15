//
// Created by lucka on 2021-07-15.
//

#ifndef LETHE_LETHEGRIDTOOLS_H
#define LETHE_LETHEGRIDTOOLS_H


#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

using namespace dealii;
namespace LetheGridTools
{
    /**
    * @brief
    * Map the vertex index to the cell that includes that vertex.
    * This map is used to find all the cell close to a specific vertex.
    */
    template <int dim>
    void
    vertices_cell_mapping(const DoFHandler<dim> &dof_handler,
            std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> & vertices_cell_map);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                     Point<dim>             &point);

    /**
      * @brief
      * Return the cell around a point based on a initial guess of a closed cell
      * (look in the neighbors of this cell)
      *
      * @param cell , The initial cell. We suspect the point of being in one of the neighbours of this cell.
      *
      * @param point, The point that we want to find the cell that contains it
      */

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree_with_guess(const DoFHandler<dim> &dof_handler,
            const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     Point<dim>             &point);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_neighbors(const DoFHandler<dim> &dof_handler,
            std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
            const typename DoFHandler<dim>::active_cell_iterator &cell,
            Point<dim>                                           &point);

    /**
      * @brief
      *Return a vector of cells around a cell including vertex neighbors
      *
      * @param cell , The initial cell. we want to know all the cells that share a vertex with this cell.
      */
    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_cell(const DoFHandler<dim> &dof_handler,
            std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
            const typename DoFHandler<dim>::active_cell_iterator &cell);



    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_around_edge(const DoFHandler<dim> &dof_handler,
                            const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                            Point<dim> &point_2);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    find_cells_in_cells(const DoFHandler<dim> &dof_handler_1,const DoFHandler<dim> &dof_handler_2,
            const typename DoFHandler<dim>::active_cell_iterator &cell);








}

#endif //LETHE_LETHEGRIDTOOLS_H

