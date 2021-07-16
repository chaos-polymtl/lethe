//
// Created by lucka on 2021-07-15.
//

#include <core/lethegridtools.h>



template <int dim>
void
LetheGridTools::vertices_cell_mapping(const DoFHandler<dim> &dof_handler,
        std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> & vertices_cell_map)
{
    // Find all the cells around each vertices

    vertices_cell_map.clear();
    const auto &cell_iterator = dof_handler.active_cell_iterators();


    // // Loop on all the cells and find their vertices to fill the map of sets of
    // cells around each vertex
    for (const auto &cell : cell_iterator)
    {
        if (cell->is_locally_owned() || cell->is_ghost())
        {
            const unsigned int vertices_per_cell =
                    GeometryInfo<dim>::vertices_per_cell;
            for (unsigned int i = 0; i < vertices_per_cell; i++)
            {
                // First obtain vertex index
                unsigned int v_index = cell->vertex_index(i);

                // Insert the cell into the set of cell around that vertex.
                vertices_cell_map[v_index].insert(cell);
            }
        }
    }
}


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

template <int dim>
typename DoFHandler<dim>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(const DoFHandler<dim> &dof_handler,
                                                      std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        Point<dim>                                            &point)
{
    // Find the cell around a point based on an initial cell.

    // Find the cells around the initial cell ( cells that share a vertex with the
    // original cell).
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
            active_neighbors_set = LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map,cell);
    // Loop over that group of cells
    for (unsigned int i = 0; i < active_neighbors_set.size(); ++i)
    {
        bool inside_cell = active_neighbors_set[i]->point_inside(point);
        if (inside_cell)
        {
            return active_neighbors_set[i];
        }
    }
    // The cell is not found near the initial cell so we use the cell tree
    // algorithm instead (much slower).
    std::cout << "Cell not found around " << point << std::endl;
    return LetheGridTools::find_cell_around_point_with_tree(dof_handler, point);
}

template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_around_cell(std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                                                        const typename DoFHandler<dim>::active_cell_iterator &cell)
{
    // Find all the cells that share a vertex with a reference cell including the
    // initial cell.
    // std::set<typename DoFHandler<dim>::active_cell_iterator> neighbors_cells;

    std::unordered_set<typename DoFHandler<dim>::active_cell_iterator, LetheGridTools::hash_cell<dim>, LetheGridTools::equal_cell<dim>> neighbors_cells;

    // Loop over the vertices of the initial cell and find all the cells around
    // each vertex and add them to the set of cells around the reference cell.
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++)
    {
        unsigned int v_index = cell->vertex_index(i);
        neighbors_cells.insert(vertices_cell_map[v_index].begin(),
                               vertices_cell_map[v_index].end());
    }
    // Transform the set into a vector.
    std::vector<typename DoFHandler<dim>::active_cell_iterator>
            cells_sharing_vertices(neighbors_cells.begin(), neighbors_cells.end());
    return cells_sharing_vertices;
}

template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_in_cells(const DoFHandler<dim> &dof_handler_1,
                    const typename DoFHandler<dim>::active_cell_iterator &cell){
    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_inside;
    for (const auto &cell_iter : dof_handler_1.active_cell_iterators())
    {
        for(unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i){
            if(cell->point_inside(cell_iter->vertex(i))){
                cells_inside.push_back(cell_iter);
                break;
            }
        }
    }

    return cells_inside;
}


template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
find_cells_around_edge(const DoFHandler<dim> &dof_handler,
                       std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                       const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                       Point<dim> &point_2){

    auto& starting_cell=find_cell_around_point_with_tree(dof_handler,point_1);
    Tensor<1,dim> unit_direction=(point_2-point_1)/(point_2-point_1).norm();
    double total_dist=(point_2-point_1).norm();

    std::set<typename DoFHandler<dim>::active_cell_iterator> cells_pierced_set;
    cells_pierced_set.push_back(starting_cell);
    Point<dim> next_point=point_1;
    auto& next_cell=starting_cell;
    double dist_done=0;
    while(dist_done<total_dist){
        next_point= next_point+unit_direction*next_cell->measure()/2;
        dist_done+=next_cell->measure()/2;
        next_cell=LetheGridTools::find_cell_around_point_with_neighbors(vertices_cell_map,next_point);
        cells_pierced_set.insert(next_cell);
    }

    std::vector<typename DoFHandler<dim>::active_cell_iterator>
            cells_pierced(cells_pierced_set.begin(), cells_pierced_set.end());
    return cells_pierced;


}

template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
find_cells_around_flat_cell(const DoFHandler<dim> &dof_handler,
                       std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                       const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                       Point<dim> &point_2){

    auto& starting_cell=find_cell_around_point_with_tree(dof_handler,point_1);
    Tensor<1,dim> unit_direction=(point_2-point_1)/(point_2-point_1).norm();
    double total_dist=(point_2-point_1).norm();

    bool all_cell_found=false;
    std::set<typename DoFHandler<dim>::active_cell_iterator> cells_pierced_set;
    cells_pierced_set.push_back(starting_cell);
    Point<dim> next_point=point_1;
    auto& next_cell=starting_cell;
    double dist_done=0;
    while(dist_done<dist_done){
        next_point= next_point+unit_direction*next_cell->measure()/2;
        dist_done+=next_cell->measure()/2;
        next_cell=LetheGridTools::find_cell_around_point_with_neighbors(vertices_cell_map,next_point);
        cells_pierced_set.insert(next_cell);
    }

    std::vector<typename DoFHandler<dim>::active_cell_iterator>
            cells_pierced(cells_pierced_set.begin(), cells_pierced_set.end());
    return cells_pierced;
}




template typename DoFHandler<3>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(const DoFHandler<3> &dof_handler,
                                 Point<3>             &point);
template typename DoFHandler<2>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(const DoFHandler<2> &dof_handler,
                                 Point<2>             &point);

template void
LetheGridTools::vertices_cell_mapping(const DoFHandler<2> &dof_handler,
                                      std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>> &vertices_cell_map);
template void
LetheGridTools::vertices_cell_mapping(const DoFHandler<3> &dof_handler,
                                      std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>> &vertices_cell_map);


template typename DoFHandler<2>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(const DoFHandler<2> &dof_handler,
                                                      std::map<unsigned int,std::set<typename DoFHandler<2>::active_cell_iterator>> &vertices_cell_map,
                                                      const typename DoFHandler<2>::active_cell_iterator &cell,
                                                      Point<2>                                            &point);

template typename DoFHandler<3>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(const DoFHandler<3> &dof_handler,
                                                      std::map<unsigned int,std::set<typename DoFHandler<3>::active_cell_iterator>> &vertices_cell_map,
                                                      const typename DoFHandler<3>::active_cell_iterator &cell,
                                                      Point<3>                                            &point);


template typename std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_cells_around_cell<2>(std::map<unsigned int,std::set<typename DoFHandler<2>::active_cell_iterator>> &vertices_cell_map,
                                       const typename DoFHandler<2>::active_cell_iterator &cell);

template typename std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_cells_around_cell<3>(std::map<unsigned int,std::set<typename DoFHandler<3>::active_cell_iterator>> &vertices_cell_map,
                                       const typename DoFHandler<3>::active_cell_iterator &cell);

template std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_cells_in_cells(const DoFHandler<2> &dof_handler_1,
                                    const typename DoFHandler<2>::active_cell_iterator &cell);

template std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_cells_in_cells(const DoFHandler<3> &dof_handler_1,
                                    const typename DoFHandler<3>::active_cell_iterator &cell);
