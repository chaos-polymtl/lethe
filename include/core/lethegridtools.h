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

#include <unordered_set>

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
    cell_pierced_by_edge(const typename DoFHandler<dim>::active_cell_iterator &cell,const typename DoFHandler<1,dim>::active_cell_iterator &cell_edge);

    template <int dim>
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
    move_grid(Triangulation<dim> mesh,Tensor<2,dim+1> displacement);


    template <int dim>
    struct hash_cell
    {
      std::size_t
      operator()(const typename DoFHandler<dim>::active_cell_iterator &cell)
        const noexcept
      {
        size_t value = 0;
        for (int i = GeometryInfo<dim>::vertices_per_cell; i > 0; --i)
          {
            value = value ^ std::hash<int>()(cell->vertex_index(i));
          }
        return value;
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
        bool is_equal = true;
        for (int i = GeometryInfo<dim>::vertices_per_cell; i > 0; --i)
          {
            if (cell1->vertex_index(i) != cell2->vertex_index(i))
              {
                is_equal = false;
                break;
              }
          }
        return is_equal;
      }
    };

}


#endif //LETHE_LETHEGRIDTOOLS_H

