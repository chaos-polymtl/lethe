//
// Created by lucka on 2021-07-15.
//

#include <core/lethegridtools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/fe/fe_q.h>
#include <cmath>


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
    //std::unordered_set<typename DoFHandler<dim>::active_cell_iterator, LetheGridTools::hash_cell<dim>, LetheGridTools::equal_cell<dim>> neighbors_cells;
    std::set<typename DoFHandler<dim>::active_cell_iterator> neighbors_cells;

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
    //Loop over all the cell of the dof handler_1 and check one of the vertices of each cell_iter is contained in cell. This loop also check if the cell_iter contained a vertex of cell. If either of those are true, the cell_iter is returned as being contained by cell.
    // The variable cell is a cell from another mesh then the mesh described by dof_handler_1.
    for (const auto &cell_iter : dof_handler_1.active_cell_iterators())
    {
        for(unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i){
            if(cell->point_inside(cell_iter->vertex(i))){
                cells_inside.push_back(cell_iter);
                break;
            }
        }
        for(unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i){
            if(cell_iter->point_inside(cell->vertex(i))){
                cells_inside.push_back(cell_iter);
                break;
            }
        }

    }

    return cells_inside;
}

template <int dim>
bool
LetheGridTools::cell_pierced_by_edge(const typename DoFHandler<dim>::active_cell_iterator &cell,const typename DoFHandler<1,dim>::active_cell_iterator &cell_edge){

    auto &local_manifold = cell_edge->get_manifold();
    std::vector<Point<dim>> manifold_points(GeometryInfo<dim-1>::vertices_per_cell);

    for (unsigned int i = 0; i < GeometryInfo<dim-1>::vertices_per_cell; ++i) {
        manifold_points[i] = cell_edge->vertex(i);
    }
    auto surrounding_points =
            make_array_view(manifold_points.begin(), manifold_points.end());
    using numbers::PI;

    // A cell that is pierced either has to:
    // A) Fill these two conditions
    //    A1) The projection of one of the cell's
    //         vertices must fall on the edge
    //    A2) At least one of the scalar product of the normal as to be negative

    bool condition_a1=false;
    bool condition_a2=false;

    for (const auto face : cell->face_indices())
    {
        auto local_face = cell->face(face);
        std::vector<Tensor<1,dim>> normals_of_face_vertex(GeometryInfo<dim>::vertices_per_face);
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i) {
            Point<dim> projected_point =
                    local_manifold.project_to_manifold(surrounding_points, local_face->vertex(i));
            if (cell_edge->point_inside(projected_point))
                condition_a1 = true;
            if(( local_face->vertex(i) - projected_point).norm()!=0)
                normals_of_face_vertex[i] =( local_face->vertex(i) - projected_point)/( local_face->vertex(i) - projected_point).norm();
            else
                return true;
        }

        for (unsigned int i = 0; i < 3; ++i) {
            if(i==2){
                Tensor<1,dim> temp =  normals_of_face_vertex[1];
                normals_of_face_vertex[1]=normals_of_face_vertex[0];
                normals_of_face_vertex[0]=temp;

            }
            if(i==3){
                Tensor<1,dim> temp =  normals_of_face_vertex[3];
                normals_of_face_vertex[3]=normals_of_face_vertex[1];
                normals_of_face_vertex[1]=temp;
            }

            for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j) {
                double s=0;
                for (unsigned int k = 0; k < 3; ++k) {
                    unsigned int index_1=j+k;
                    unsigned int index_2=j+k+1;
                    if((j+k)>=GeometryInfo<dim>::vertices_per_face){
                        index_1=j+k-GeometryInfo<dim>::vertices_per_face;
                    }
                    if((j+k+1)>=GeometryInfo<dim>::vertices_per_face){
                        index_1=j+k+1-GeometryInfo<dim>::vertices_per_face;
                    }

                    s+=std::acos( scalar_product(normals_of_face_vertex[index_1], normals_of_face_vertex[index_2])/(normals_of_face_vertex[index_1].norm()*normals_of_face_vertex[index_2].norm()));
                }
                if(s<PI){
                    condition_a2 = true;
                }
            }
        }
        if (condition_a1 && condition_a2)
            return true;

    }



    return false;
}


template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_around_edge(const DoFHandler<dim> &dof_handler,
                       std::map<unsigned int,std::set<typename DoFHandler<dim>::active_cell_iterator>> &vertices_cell_map,
                       const typename DoFHandler<dim>::active_cell_iterator &cell,Point<dim> &point_1,
                       Point<dim> &point_2){




    std::set<typename DoFHandler<dim>::active_cell_iterator> cells_pierced_set;
    DoFHandler<1, dim>    local_edge_dof_handler;
    Triangulation< 1, dim> local_edge_triangulation;
    std::vector<Point<dim>>        vertices_of_edge(2);
    std::vector<CellData<1>> local_edge_cell_data(1);
    FESystem<dim - 1, dim> local_face_fe(
            FE_Q<dim - 1, dim>(1));

    vertices_of_edge[0]=point_1;
    vertices_of_edge[1]=point_2;
    local_edge_cell_data[0].vertices[0] = 0;
    local_edge_cell_data[0].vertices[1] = 1;

    local_edge_triangulation.create_triangulation(
            vertices_of_edge,
            local_edge_cell_data,
            SubCellData());
    local_edge_dof_handler.reinit(
            local_edge_triangulation);
    local_edge_dof_handler.distribute_dofs(local_face_fe);

    const typename DoFHandler<1,dim>::active_cell_iterator edge_cell;
    edge_cell=local_edge_dof_handler.active_cell_iterators().begin();
    if(dim==2){

        cells_pierced_set = find_cells_around_flat_cell(dof_handler, edge_cell, vertices_cell_map);

    }
    if(dim==3){
        auto& starting_cell=find_cell_around_point_with_tree(dof_handler,point_1);
        Tensor<1,dim> unit_direction=(point_2-point_1)/(point_2-point_1).norm();
        double total_dist=(point_2-point_1).norm();

        std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
                LetheGridTools::hash_cell<dim>,
                LetheGridTools::equal_cell<dim>>
                previous_candidate_cells;

        std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
                LetheGridTools::hash_cell<dim>,
                LetheGridTools::equal_cell<dim>>
                current_candidate_cells;

        std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
                LetheGridTools::hash_cell<dim>,
                LetheGridTools::equal_cell<dim>>
                intersected_cells;

        previous_candidate_cells.insert(starting_cell);
        current_candidate_cells.insert(starting_cell);
        intersected_cells.insert(starting_cell);

        int n_previous_intersected = 0;

        while (intersected_cells.size() > n_previous_intersected) {
            n_previous_intersected=intersected_cells.size();
            // Find all cells around previous candidate cells
            for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
                    previous_candidate_cells) {
                current_candidate_cells.insert(
                        LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map,
                                                                    cell_iter));
            }

            // Reset the list of previous candidates
            previous_candidate_cells.clear();

            // Check if current candidate cells are intersected, and if they are, add
            // them to the intersected cell set. If the added cell was not already
            // present in the intersected cell set, add it to the previous candidate
            // cells as well.
            for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
                    current_candidate_cells) {
                if (LetheGridTools::cell_pierced_by_edge(cell_iter, edge_cell))
                {
                    // If the cell was not present in the intersected cells set
                    if (intersected_cells.insert(cell_iter).second) {
                        previous_candidate_cells.insert(cell_iter);
                    }
                }
            }
            current_candidate_cells.clear();
        }





    }


    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_pierced(cells_pierced_set.begin(),cells_pierced_set.end());

    return cells_pierced;
}

template <int dim>
bool
LetheGridTools::cell_cut_by_flat(
        const typename DoFHandler<dim>::active_cell_iterator &cell,
        const typename DoFHandler<dim-1, dim>::active_cell_iterator &cell_flat) {
    auto &local_manifold = cell_flat->get_manifold();
    std::vector<Point<dim>> manifold_points(GeometryInfo<dim-1>::vertices_per_cell);

    for (unsigned int i = 0; i < GeometryInfo<dim-1>::vertices_per_cell; ++i) {
        manifold_points[i] = cell_flat->vertex(i);
    }
    auto surrounding_points =
            make_array_view(manifold_points.begin(), manifold_points.end());

    // A cell that is cut either has to:
    // A) Contain a vertex from the flat
    // B) Fill these two conditions
    //    B1) The projection of one of the cell's
    //         vertices must fall on the flat
    //    B2) The cell must have vertices on each side of the flat
    // C) Fill these conditions
    //    C1) One of the faces of the cell must have vertices of the flat cell on either side of it.
    //    C2) The cell must have vertices on each side of the flat (identical to B2)


    // Check for condition A
    for (const Point<dim> &point : manifold_points) {
        if (cell->point_inside(point))
            return true;
    }

    // Check for condition B
    bool condition_B1 = false;
    bool condition_B2 = false;
    Tensor<1, dim> last_normal;
    last_normal.clear();

    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i) {
        Point<dim> projected_point =
                local_manifold.project_to_manifold(surrounding_points, cell->vertex(i));
        Tensor<1, dim> normal = cell->vertex(i) - projected_point;

        // Check if the projected vertex falls inside the flat
        if (cell->point_inside(projected_point))
            condition_B1 = true;

        // Check if we switched to the other side
        // of the flat during this iteration
        double scalar_prod = scalar_product(normal, last_normal);
        last_normal = normal;
        if (scalar_prod < 0)
            condition_B2 = true;

        if (condition_B1 && condition_B2)
            return true;
    }


    if(condition_B2) {
        bool condition_C2 = false;
        manifold_points.clear();
        for (const auto face : cell->face_indices()) {
            auto face_iter = cell->face(face);
            last_normal.clear();
            auto &local_face_manifold = face_iter->get_manifold();

            for (unsigned int i = 0; i < GeometryInfo<dim-1>::vertices_per_cell; ++i) {
                manifold_points[i] = face_iter->vertex(i);
            }
            auto surrounding_points_face =
                    make_array_view(manifold_points.begin(), manifold_points.end());

            for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i) {
                Point<dim> projected_point =
                        local_face_manifold.project_to_manifold(surrounding_points_face, cell_flat->vertex(i));
                Tensor<1, dim> normal = cell->vertex(i) - projected_point;
                double scalar_prod = scalar_product(normal, last_normal);
                last_normal = normal;
                if (scalar_prod < 0)
                    condition_C2 = true;

            }
            if (condition_B1 && condition_C2)
                return true;
        }
    }
    return false;
}


template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_around_flat_cell(
  const DoFHandler<dim> &                                   dof_handler,
  const typename DoFHandler<dim - 1,dim>::active_cell_iterator &cell,
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &vertices_cell_map)
{
    std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_cut;
      auto &starting_cell =
              find_cell_around_point_with_tree(dof_handler, cell->vertex(0));

      std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
              LetheGridTools::hash_cell<dim>,
              LetheGridTools::equal_cell<dim>>
              previous_candidate_cells;

      std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
              LetheGridTools::hash_cell<dim>,
              LetheGridTools::equal_cell<dim>>
              current_candidate_cells;

      std::unordered_set<typename DoFHandler<dim>::active_cell_iterator,
              LetheGridTools::hash_cell<dim>,
              LetheGridTools::equal_cell<dim>>
              intersected_cells;

      previous_candidate_cells.insert(starting_cell);
      current_candidate_cells.insert(starting_cell);
      intersected_cells.insert(starting_cell);

      int n_previous_intersected = 0;

      while (intersected_cells.size() > n_previous_intersected) {
          n_previous_intersected=intersected_cells.size();
          // Find all cells around previous candidate cells
          for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
                  previous_candidate_cells) {
              current_candidate_cells.insert(
                      LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map,
                                                                  cell_iter));
          }

          // Reset the list of previous candidates
          previous_candidate_cells.clear();

          // Check if current candidate cells are intersected, and if they are, add
          // them to the intersected cell set. If the added cell was not already
          // present in the intersected cell set, add it to the previous candidate
          // cells as well.
          for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
                  current_candidate_cells) {
              if (LetheGridTools::cell_cut_by_flat(cell_iter, cell))
                {
                  // If the cell was not present in the intersected cells set
                  if (intersected_cells.insert(cell_iter).second) {
                      previous_candidate_cells.insert(cell_iter);
                  }
              }
          }
          current_candidate_cells.clear();
      }
      cells_cut(intersected_cells.begin(), intersected_cells.end());


  return cells_cut;
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

template bool
LetheGridTools::cell_cut_by_flat<2>(
        const typename DoFHandler<2>::active_cell_iterator &cell,
        const typename DoFHandler<2-1, 2>::active_cell_iterator &cell_flat);

template bool
LetheGridTools::cell_cut_by_flat<3>(
        const typename DoFHandler<3>::active_cell_iterator &cell,
        const typename DoFHandler<3-1, 3>::active_cell_iterator &cell_flat);
