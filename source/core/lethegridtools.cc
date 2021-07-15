//
// Created by lucka on 2021-07-15.
//

#include <core/lethegridtools.h>


template <int dim>
typename DoFHandler<dim>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                 Point<dim>             &point)
{
    // find cell around point using tree properties instead of looping on all cell

    // define stuff for the search
    MappingQ1<dim> map;
    const auto &   cell_iterator = dof_handler.cell_iterators_on_level(0);
    unsigned int   max_childs    = GeometryInfo<dim>::max_children_per_cell;
    typename DoFHandler<dim>::cell_iterator best_cell_iter;

    bool cell_0_found = false;

    // loop on the cell on lvl 0 of the mesh
    for (const auto &cell : cell_iterator)
    {
        try
        {
            const Point<dim, double> p_cell =
                    map.transform_real_to_unit_cell(cell, point);

            const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

            if (dist == 0)
            {
                // cell on lvl 0 found
                cell_0_found   = true;
                best_cell_iter = cell;
            }
        }
        catch (const typename MappingQGeneric<dim>::ExcTransformationFailed &)
        {}
    }

    if (cell_0_found)
    {
        // the cell on lvl 0 contain the point so now loop on the childs of
        // this cell when we found the child of the cell that containt it we
        // stop and loop if the cell is active if the cell is not active loop
        // on the child of the cell. Repeat
        unsigned int lvl = 0;
        while (best_cell_iter->is_active() == false)
        {
            bool         cell_found = false;
            double       best_dist  = DBL_MAX;
            unsigned int best_index = 0;
            for (unsigned int i = 0; i < max_childs; ++i)
            {
                try
                {
                    const Point<dim, double> p_cell =
                            map.transform_real_to_unit_cell(best_cell_iter->child(i),
                                                            point);
                    const double dist =
                            GeometryInfo<dim>::distance_to_unit_cell(p_cell);
                    bool inside = true;


                    if (dist <= best_dist and inside)
                    {
                        best_dist  = dist;
                        best_index = i;
                        cell_found = true;
                        if (dist == 0)
                            break;
                    }
                }
                catch (
                        const typename MappingQGeneric<dim>::ExcTransformationFailed &)
                {}
            }

            best_cell_iter = best_cell_iter->child(best_index);

            if (cell_found == false)
            {
                std::cout << "cell not found" << point << std::endl;
                break;
            }

            lvl += 1;
        }
    }


    return best_cell_iter;
}

template typename DoFHandler<3>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(const DoFHandler<3> &dof_handler,
                                 Point<3>             &point);
template typename DoFHandler<2>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(const DoFHandler<2> &dof_handler,
                                 Point<2>             &point);
