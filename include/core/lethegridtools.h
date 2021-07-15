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
    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                     Point<dim>             &point);


    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_tree_with_guess(const DoFHandler<dim> &dof_handler, const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     Point<dim>             &point);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterator
    find_cell_around_point_with_neighbors(const DoFHandler<dim> &dof_handler,
            const typename DoFHandler<dim>::active_cell_iterator &cell,
            Point<dim>                                            &point);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterators
    find_cells_around_cells(const DoFHandler<dim> &dof_handler,
            const typename DoFHandler<dim>::active_cell_iterator &cell);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterators
    find_cells_around_edge(const DoFHandler<dim> &dof_handler,
                            const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                            Point<dim> &point_2);

    template <int dim>
    typename DoFHandler<dim>::active_cell_iterators
    find_cells_in_cells(const DoFHandler<dim> &dof_handler_1,const DoFHandler<dim> &dof_handler_2,
            const typename DoFHandler<dim>::active_cell_iterator &cell);








}

#endif //LETHE_LETHEGRIDTOOLS_H

