#include <core/lethegridtools.h>
#include <core/serial_solid.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/dem_properties.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#include <cmath>
#include <unordered_map>


template <int dim>
void
LetheGridTools::vertices_cell_mapping(
  const DoFHandler<dim> &dof_handler,
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &vertices_cell_map)
{
  vertices_cell_map.clear();
  const auto &cell_iterator = dof_handler.active_cell_iterators();

  // Loop on all the cells and find their vertices, and use them to fill
  // the map of sets of cells around each vertex
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
LetheGridTools::find_cell_around_point_with_tree(
  const DoFHandler<dim> &dof_handler,
  const Point<dim> &     point)
{
  // Define temporary variables to store search parameters and intermediate
  // results
  MappingQ1<dim> mapping;
  const auto &   cell_iterator = dof_handler.cell_iterators_on_level(0);

  typename DoFHandler<dim>::cell_iterator best_cell_iter;

  bool cell_on_level_0_found = false;

  // Loop on the cells on level 0 of the mesh
  for (const auto &cell : cell_iterator)
    {
      try
        {
          const Point<dim, double> p_cell =
            mapping.transform_real_to_unit_cell(cell, point);

          const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

          if (dist < 1e-12)
            {
              // cell on lvl 0 found
              cell_on_level_0_found = true;
              best_cell_iter        = cell;
              break;
            }
        }
      catch (const typename MappingQGeneric<dim>::ExcTransformationFailed &)
        {}
    }
  double best_dist_last = DBL_MAX;
  if (cell_on_level_0_found)
    {
      // A cell on level 0 contains the point, so we loop on the children of
      // this cell. When we found the child of the cell that contains the point,
      // we check if the cell is active and stop if it is. Otherwise, we repeat
      // this process until we find the active cell that contains the point.

      unsigned int lvl        = 0;
      unsigned int max_childs = GeometryInfo<dim>::max_children_per_cell;
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
                    mapping.transform_real_to_unit_cell(
                      best_cell_iter->child(i), point);
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
              break;
            }

          lvl += 1;
          best_dist_last = best_dist;
        }
    }

  if (best_dist_last >= 1e-9 && cell_on_level_0_found == false)
    {
      throw std::runtime_error("The point is not inside the mesh");
    }
  return best_cell_iter;
}


template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_boundary_cells_in_sphere(
  const DoFHandler<dim> &dof_handler,
  const Point<dim> &     center,
  const double           radius)
{
  MappingQ1<dim> mapping;
  const auto &   cell_iterator = dof_handler.cell_iterators_on_level(0);
  unsigned int   max_childs    = GeometryInfo<dim>::max_children_per_cell;

  std::set<typename DoFHandler<dim>::cell_iterator> boundary_cells_candidates;
  std::set<typename DoFHandler<dim>::cell_iterator>
                                                    last_boundary_cells_candidates;
  std::set<typename DoFHandler<dim>::cell_iterator> cells_at_boundary;
  // Loop over the cells on level 0 of the mesh
  for (const auto &cell : cell_iterator)
    {
      if (cell->at_boundary())
        {
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell;
               ++i)
            {
              if ((cell->vertex(i) - center).norm() <= radius ||
                  cell->point_inside(center))
                {
                  boundary_cells_candidates.insert(cell);
                  break;
                }
            }
        }
    }
  // Loop over the boundary cells on level 0 to check their children.
  last_boundary_cells_candidates = boundary_cells_candidates;
  if (boundary_cells_candidates.size() != 0)
    {
      bool all_cell_are_active = false;
      while (!all_cell_are_active)
        {
          // Only stop when all the candidate cell are active
          all_cell_are_active = true;
          boundary_cells_candidates.clear();
          // loop over the last set of candidates.
          for (auto cell : last_boundary_cells_candidates)
            {
              if (cell->at_boundary())
                {
                  for (unsigned int i = 0;
                       i < GeometryInfo<dim>::vertices_per_cell;
                       ++i)
                    {
                      // Check if the cell is in the sphere.
                      if ((cell->vertex(i) - center).norm() <= radius ||
                          cell->point_inside(center))
                        {
                          if (cell->is_active())
                            {
                              cells_at_boundary.insert(cell);
                            }
                          else
                            {
                              // If we are here, the cell has children.
                              all_cell_are_active = false;
                              for (unsigned int j = 0; j < max_childs; ++j)
                                {
                                  // Store the children of the cell for further
                                  // checks.
                                  boundary_cells_candidates.insert(
                                    cell->child(j));
                                }
                            }
                          break;
                        }
                    }
                }
            }
          last_boundary_cells_candidates.clear();
          last_boundary_cells_candidates = boundary_cells_candidates;
        }
    }
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    cells_at_boundary_v(cells_at_boundary.begin(), cells_at_boundary.end());
  // Output all the active cells in the sphere radius that are at a boundary.
  return cells_at_boundary_v;
}


template <int dim>
typename DoFHandler<dim>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(
  const DoFHandler<dim> &dof_handler,
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &                                                   vertices_cell_map,
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const Point<dim> &                                    point)
{
  // Find the cells that share a vertex with the original cell.
  std::vector<typename DoFHandler<dim>::active_cell_iterator>
    active_neighbors_set =
      LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map, cell);
  // Loop over that group of cells
  for (unsigned int i = 0; i < active_neighbors_set.size(); ++i)
    {
      bool inside_cell = active_neighbors_set[i]->point_inside(point);
      if (inside_cell)
        {
          return active_neighbors_set[i];
        }
    }
  // The cell is not found near the initial cell, so we use the cell tree
  // algorithm instead (much slower).
  return LetheGridTools::find_cell_around_point_with_tree(dof_handler, point);
}

template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_around_cell(
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &                                                   vertices_cell_map,
  const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  // Find all the cells that share a vertex with a reference cell, including the
  // initial cell.
  std::set<typename DoFHandler<dim>::active_cell_iterator> neighbors_cells;

  // Loop over the vertices of the initial cell, find all the cells around
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
LetheGridTools::find_cells_in_cells(
  const DoFHandler<dim> &                               dof_handler,
  const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_inside;
  // Loop over all the cells of the dof handler and check one of the vertices
  // of each cell is contained in the external cell, or vice versa. If either
  // of those are true, the cell_iter is added to the set of cells contained by
  // the external cell.

  for (const auto &cell_iter : dof_handler.active_cell_iterators())
    {
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          if (cell->point_inside(cell_iter->vertex(i)))
            {
              cells_inside.push_back(cell_iter);
              break;
            }
        }
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          if (cell_iter->point_inside(cell->vertex(i)))
            {
              cells_inside.push_back(cell_iter);
              break;
            }
        }
    }

  return cells_inside;
}


template <int dim>
bool
LetheGridTools::cell_pierced_by_edge(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const TriaIterator<CellAccessor<1, dim>> &            cell_edge)
{
  std::vector<Point<dim>> manifold_points(GeometryInfo<1>::vertices_per_cell);

  for (unsigned int i = 0; i < GeometryInfo<1>::vertices_per_cell; ++i)
    {
      manifold_points[i] = cell_edge->vertex(i);
    }

  using numbers::PI;

  // A cell that is pierced either has to fill these two conditions:
  //    A1) The projection of one of the cell's
  //         vertices must fall on the edge
  //    A2) At least one of the scalar product of the normal as to be negative

  for (const auto face : cell->face_indices())
    {
      bool                        condition_a1 = false;
      bool                        condition_a2 = true;
      auto                        local_face   = cell->face(face);
      std::vector<Tensor<1, dim>> normals_of_face_vertex(
        GeometryInfo<dim>::vertices_per_face);
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
        {
          Point<dim> projected_point =
            GridTools::project_to_object(cell_edge, local_face->vertex(i));

          if (cell_edge->point_inside(projected_point) &&
              cell->point_inside(projected_point))
            {
              condition_a1 = true;
            }
          if ((local_face->vertex(i) - projected_point).norm() != 0)
            normals_of_face_vertex[i] =
              (local_face->vertex(i) - projected_point) /
              (local_face->vertex(i) - projected_point).norm();
          else
            return true;
        }

      for (unsigned int i = 0; i < 3; ++i)
        {
          if (i == 2)
            {
              Tensor<1, dim> temp       = normals_of_face_vertex[1];
              normals_of_face_vertex[1] = normals_of_face_vertex[0];
              normals_of_face_vertex[0] = temp;
            }
          if (i == 3)
            {
              Tensor<1, dim> temp       = normals_of_face_vertex[3];
              normals_of_face_vertex[3] = normals_of_face_vertex[1];
              normals_of_face_vertex[1] = temp;
            }

          for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face;
               ++j)
            {
              double s = 0;
              for (unsigned int k = 0; k < 3; ++k)
                {
                  unsigned int index_1 = j + k;
                  unsigned int index_2 = j + k + 1;
                  if ((j + k) >= GeometryInfo<dim>::vertices_per_face)
                    {
                      index_1 = j + k - GeometryInfo<dim>::vertices_per_face;
                    }
                  if ((j + k + 1) >= GeometryInfo<dim>::vertices_per_face)
                    {
                      index_2 =
                        j + k + 1 - GeometryInfo<dim>::vertices_per_face;
                    }

                  double dot = scalar_product(normals_of_face_vertex[index_1],
                                              normals_of_face_vertex[index_2]);
                  s += std::acos(std::clamp(dot, -1.0, 1.0));
                }
              if (s < PI)
                {
                  condition_a2 = false;
                  break;
                }
            }
        }
      if (condition_a2)
        if (condition_a1 && condition_a2)
          {
            return true;
          }
    }



  return false;
}

template <int dim>
bool
LetheGridTools::cell_pierced_by_edge(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  Point<dim>                                            point_1,
  Point<dim>                                            point_2)
{
  Triangulation<1, dim>    local_edge_triangulation;
  std::vector<Point<dim>>  vertices_of_edge(2);
  std::vector<CellData<1>> local_edge_cell_data(1);

  vertices_of_edge[0]                 = point_1;
  vertices_of_edge[1]                 = point_2;
  local_edge_cell_data[0].vertices[0] = 0;
  local_edge_cell_data[0].vertices[1] = 1;

  local_edge_triangulation.create_triangulation(vertices_of_edge,
                                                local_edge_cell_data,
                                                SubCellData());

  for (const auto &edge_cell : local_edge_triangulation.active_cell_iterators())
    {
      return LetheGridTools::cell_pierced_by_edge<dim>(cell, edge_cell);
    }

  // Default return statement in case loop is empty;
  return false;
}



template <int spacedim, int structdim>
std::pair<std::pair<Point<spacedim>, bool>, Tensor<1, spacedim>>
LetheGridTools::project_to_d_linear_object(
  const typename DoFHandler<structdim, spacedim>::active_cell_iterator &object,
  const Point<spacedim> &trial_point)
{
  // let's look at this for simplicity for a quad (structdim==2) in a
  // space with spacedim>2 (notate trial_point by y): all points on the
  // surface are given by
  //   x(\xi) = sum_i v_i phi_x(\xi)
  // where v_i are the vertices of the quad, and \xi=(\xi_1,\xi_2) are the
  // reference coordinates of the quad. so what we are trying to do is
  // find a point x on the surface that is closest to the point y. there
  // are different ways to solve this problem, but in the end it's a
  // nonlinear problem and we have to find reference coordinates \xi so
  // that J(\xi) = 1/2 || x(\xi)-y ||^2 is minimal. x(\xi) is a function
  // that is structdim-linear in \xi, so J(\xi) is a polynomial of degree
  // 2*structdim that we'd like to minimize. unless structdim==1, we'll
  // have to use a Newton method to find the answer. This leads to the
  // following formulation of Newton steps:
  //
  // Given \xi_k, find \delta\xi_k so that
  //   H_k \delta\xi_k = - F_k
  // where H_k is an approximation to the second derivatives of J at
  // \xi_k, and F_k is the first derivative of J.  We'll iterate this a
  // number of times until the right hand side is small enough. As a
  // stopping criterion, we terminate if ||\delta\xi||<eps.
  //
  // As for the Hessian, the best choice would be
  //   H_k = J''(\xi_k)
  // but we'll opt for the simpler Gauss-Newton form
  //   H_k = A^T A
  // i.e.
  //   (H_k)_{nm} = \sum_{i,j} v_i*v_j *
  //                   \partial_n phi_i *
  //                   \partial_m phi_j
  // we start at xi=(0.5, 0.5).
  Point<structdim> xi;
  for (unsigned int d = 0; d < structdim; ++d)
    xi[d] = 0.5;

  Point<spacedim> x_k;
  for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
    x_k += object->vertex(i) *
           GeometryInfo<structdim>::d_linear_shape_function(xi, i);

  do
    {
      Tensor<1, structdim> F_k;
      for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
        F_k += (x_k - trial_point) * object->vertex(i) *
               GeometryInfo<structdim>::d_linear_shape_function_gradient(xi, i);

      Tensor<2, structdim> H_k;
      for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
        for (const unsigned int j : GeometryInfo<structdim>::vertex_indices())
          {
            Tensor<2, structdim> tmp = outer_product(
              GeometryInfo<structdim>::d_linear_shape_function_gradient(xi, i),
              GeometryInfo<structdim>::d_linear_shape_function_gradient(xi, j));
            H_k += (object->vertex(i) * object->vertex(j)) * tmp;
          }

      const Tensor<1, structdim> delta_xi = -invert(H_k) * F_k;
      xi += delta_xi;

      x_k = Point<spacedim>();
      for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
        x_k += object->vertex(i) *
               GeometryInfo<structdim>::d_linear_shape_function(xi, i);

      if (delta_xi.norm() < 1e-7)
        {
          break;
        }
    }
  while (true);
  Tensor<1, spacedim> normal;
  if (spacedim == 3)
    {
      double           dx    = 1e-4;
      Point<structdim> xi_dx = xi;
      Point<structdim> xi_dy = xi;
      xi_dx[0]               = xi[0] + dx;
      xi_dy[1]               = xi[1] - dx;
      Point<spacedim> x_k_dx(0, 0, 0);
      Point<spacedim> x_k_dy(0, 0, 0);
      for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
        {
          x_k_dx += object->vertex(i) *
                    GeometryInfo<structdim>::d_linear_shape_function(xi_dx, i);
          x_k_dy += object->vertex(i) *
                    GeometryInfo<structdim>::d_linear_shape_function(xi_dy, i);
        }

      normal = cross_product_3d((x_k_dx - x_k), (x_k_dy - x_k));
      normal = normal / normal.norm();
    }

  if (spacedim == 2)
    {
      Tensor<1, spacedim> temp = object->vertex(1) - object->vertex(0);
      normal[0]                = -temp[1];
      normal[1]                = temp[0];
    }
  bool inside = true;
  for (unsigned int i = 0; i < structdim; ++i)
    {
      if (xi[i] < 0 || xi[i] > 1)
        inside = false;
    }

  std::pair<Point<spacedim>, bool> point = std::make_pair(x_k, inside);
  std::pair<std::pair<Point<spacedim>, bool>, Tensor<1, spacedim>> output =
    std::make_pair(point, normal);


  return output;
}


template <int dim>
std::vector<typename DoFHandler<dim>::active_cell_iterator>
LetheGridTools::find_cells_around_edge(
  const DoFHandler<dim> &dof_handler,
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &                                                   vertices_cell_map,
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  Point<dim> &                                          point_1,
  Point<dim> &                                          point_2)
{
  std::set<typename DoFHandler<dim>::active_cell_iterator> cells_pierced_set;
  DoFHandler<1, dim>       local_edge_dof_handler;
  Triangulation<1, dim>    local_edge_triangulation;
  std::vector<Point<dim>>  vertices_of_edge(2);
  std::vector<CellData<1>> local_edge_cell_data(1);
  FESystem<dim - 1, dim>   local_face_fe(FE_Q<dim - 1, dim>(1));

  vertices_of_edge[0]                 = point_1;
  vertices_of_edge[1]                 = point_2;
  local_edge_cell_data[0].vertices[0] = 0;
  local_edge_cell_data[0].vertices[1] = 1;

  local_edge_triangulation.create_triangulation(vertices_of_edge,
                                                local_edge_cell_data,
                                                SubCellData());
  local_edge_dof_handler.reinit(local_edge_triangulation);
  local_edge_dof_handler.distribute_dofs(local_face_fe);

  const typename DoFHandler<1, dim>::active_cell_iterator edge_cell;
  edge_cell = local_edge_dof_handler.active_cell_iterators().begin();
  if (dim == 2)
    {
      cells_pierced_set =
        find_cells_around_flat_cell(dof_handler, edge_cell, vertices_cell_map);
    }
  if (dim == 3)
    {
      auto &starting_cell =
        find_cell_around_point_with_tree(dof_handler, point_1);
      Tensor<1, dim> unit_direction =
        (point_2 - point_1) / (point_2 - point_1).norm();
      double total_dist = (point_2 - point_1).norm();

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

      while (intersected_cells.size() > n_previous_intersected)
        {
          n_previous_intersected = intersected_cells.size();
          // Find all cells around previous candidate cells
          for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
               previous_candidate_cells)
            {
              current_candidate_cells.insert(
                LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map,
                                                            cell_iter));
            }

          // Reset the list of previous candidates
          previous_candidate_cells.clear();

          // Check if current candidate cells are intersected, and if they are,
          // add them to the intersected cell set. If the added cell was not
          // already present in the intersected cell set, add it to the previous
          // candidate cells as well.
          for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
               current_candidate_cells)
            {
              if (LetheGridTools::cell_pierced_by_edge(cell_iter, edge_cell))
                {
                  // If the cell was not present in the intersected cells set
                  if (intersected_cells.insert(cell_iter).second)
                    {
                      previous_candidate_cells.insert(cell_iter);
                    }
                }
            }
          current_candidate_cells.clear();
        }
    }


  std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_pierced(
    cells_pierced_set.begin(), cells_pierced_set.end());

  return cells_pierced;
}

template <int dim>
bool
LetheGridTools::cell_cut_by_flat(
  const typename DoFHandler<dim>::active_cell_iterator &         cell,
  const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell_flat)
{
  std::vector<Point<dim>> manifold_points(
    GeometryInfo<dim - 1>::vertices_per_cell);
  for (unsigned int i = 0; i < GeometryInfo<dim - 1>::vertices_per_cell; ++i)
    {
      manifold_points[i] = cell_flat->vertex(i);
    }

  // A cell that is cut either has to:
  // A) Contain a vertex from the flat
  // B) Fill these two conditions
  //    B1) The projection of one of the cell's
  //         vertices must fall on the flat
  //    B2) The cell must have vertices on each side of the flat
  // C) Fill these conditions
  //    C1) One of the faces of the cell must have vertices of the flat cell on
  //    either side of it.
  //    C2) The cell must have vertices on each side of the
  //    flat (identical to B2)


  // Check for condition A
  for (unsigned int i = 0; i < GeometryInfo<dim - 1>::vertices_per_cell; ++i)
    {
      if (cell->point_inside(cell_flat->vertex(i)))
        return true;
    }

  if constexpr (dim == 3)
    {
      // Check for condition A
      for (const auto face : cell_flat->face_indices())
        {
          auto local_face = cell_flat->face(face);
          if (LetheGridTools::cell_pierced_by_edge<3>(cell,
                                                      local_face->vertex(0),
                                                      local_face->vertex(1)))
            return true;
        }
    }

  // Check for condition B
  bool           condition_B1 = false;
  bool           condition_B2 = false;
  Tensor<1, dim> last_normal;
  last_normal.clear();
  double last_scalar_prod = 0;

  for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      auto projected_point =
        LetheGridTools::project_to_d_linear_object<dim, dim - 1>(
          cell_flat, cell->vertex(i));

      Tensor<1, dim> normal = cell->vertex(i) - projected_point.first.first;

      // Check if the projected vertex falls inside the flat

      if (projected_point.first.second &&
          cell->point_inside(projected_point.first.first))
        {
          condition_B1 = true;
        }

      // Check if we switched to the other side
      // of the flat during this iteration
      double scalar_prod = scalar_product(normal, projected_point.second);

      if (scalar_prod * last_scalar_prod < 0)
        {
          condition_B2 = true;
        }
      last_scalar_prod = scalar_prod;
      if (condition_B1 && condition_B2)
        {
          return true;
        }
    }


  if (condition_B2)
    {
      bool condition_C2 = false;
      manifold_points.clear();
      for (const auto face : cell->face_indices())
        {
          auto face_iter = cell->face(face);
          last_normal.clear();
          auto &local_face_manifold = face_iter->get_manifold();

          for (unsigned int i = 0; i < GeometryInfo<dim - 1>::vertices_per_cell;
               ++i)
            {
              manifold_points[i] = face_iter->vertex(i);
            }
          auto surrounding_points_face =
            make_array_view(manifold_points.begin(), manifold_points.end());

          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face;
               ++i)
            {
              Point<dim> projected_point =
                local_face_manifold.project_to_manifold(surrounding_points_face,
                                                        cell_flat->vertex(i));
              Tensor<1, dim> normal      = cell->vertex(i) - projected_point;
              double         scalar_prod = scalar_product(normal, last_normal);
              last_normal                = normal;
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
  const DoFHandler<dim> &                                        dof_handler,
  const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell,
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    &vertices_cell_map)
{
  TriaActiveIterator<DoFCellAccessor<dim, dim, 0>> starting_cell =
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

  size_t n_previous_intersected = 0;

  while (intersected_cells.size() > n_previous_intersected)
    {
      n_previous_intersected = intersected_cells.size();
      // Find all cells around previous candidate cells
      for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
           previous_candidate_cells)
        {
          auto new_cells =
            LetheGridTools::find_cells_around_cell<dim>(vertices_cell_map,
                                                        cell_iter);
          for (unsigned int i = 0; i < new_cells.size(); ++i)
            {
              current_candidate_cells.insert(new_cells[i]);
            }
        }

      // Reset the list of previous candidates
      previous_candidate_cells.clear();

      // Check if current candidate cells are intersected, and if they are, add
      // them to the set  of intersected cells . If the added cells were not
      // already present in the set of intersected cells, add them to the
      // previous candidate cells as well.
      for (const typename DoFHandler<dim>::active_cell_iterator &cell_iter :
           current_candidate_cells)
        {
          if (LetheGridTools::cell_cut_by_flat<dim>(cell_iter, cell))
            {
              previous_candidate_cells.insert(cell_iter);
              intersected_cells.insert(cell_iter);
              // If the cell was not present in the intersected cells set
            }
        }
      current_candidate_cells.clear();
    }

  std::vector<typename DoFHandler<dim>::active_cell_iterator> cells_cut(
    intersected_cells.begin(), intersected_cells.end());


  return cells_cut;
}

template <int spacedim, int structdim>
std::map<
  typename DoFHandler<spacedim>::active_cell_iterator,
  std::map<unsigned int,
           typename DoFHandler<structdim, spacedim>::active_cell_iterator>>
LetheGridTools::find_cells_cut_by_object(
  const DoFHandler<spacedim> &dof_handler,
  std::map<unsigned int,
           std::set<typename DoFHandler<spacedim>::active_cell_iterator>>
    &                                            vertices_cell_map,
  std::vector<SerialSolid<structdim, spacedim>> &list_of_objects)
{
  std::map<
    typename DoFHandler<spacedim>::active_cell_iterator,
    std::map<unsigned int,
             typename DoFHandler<structdim, spacedim>::active_cell_iterator>>
    cells_cut_by_object;
  if constexpr (structdim == spacedim - 1)
    {
      for (unsigned int i = 0; i < list_of_objects.size(); ++i)
        {
          const auto &object = list_of_objects[i].get_solid_dof_handler();

          const auto &object_cell_iterator = object.active_cell_iterators();
          // Loop on all the cells and find their vertices, and use them to fill
          // the map of sets of cells around each vertex
          for (const auto &cell : object_cell_iterator)
            {
              std::vector<typename DoFHandler<spacedim>::active_cell_iterator>
                cells_cut = LetheGridTools::find_cells_around_flat_cell(
                  dof_handler, cell, vertices_cell_map);
              for (unsigned int j = 0; j < cells_cut.size(); ++j)
                {
                  cells_cut_by_object[cells_cut[j]]
                                     [list_of_objects[i].get_solid_id()] = cell;
                }
            }
        }
      return cells_cut_by_object;
    }
  if constexpr (structdim == spacedim)
    {
      for (unsigned int i = 0; i < list_of_objects.size(); ++i)
        {
          auto &object = list_of_objects[i].get_solid_dof_handler();

          const auto &object_cell_iterator = object.active_cell_iterators();
          // Loop on all the cells and find their vertices, and use them to fill
          // the map of sets of cells around each vertex
          for (const auto &cell : object_cell_iterator)
            {
              if (cell->at_boundary())
                {
                  for (const auto face : cell->face_indices())
                    {
                      auto local_face = cell->face(face);
                      if (cell->at_boundary())
                        {
                          std::vector<
                            typename DoFHandler<spacedim>::active_cell_iterator>
                            cells_cut =
                              LetheGridTools::find_cells_around_flat_cell(
                                dof_handler, cell, vertices_cell_map);
                          for (unsigned int j = 0; j < cells_cut.size(); ++j)
                            {
                              cells_cut_by_object[cells_cut[j]]
                                                 [list_of_objects[i]
                                                    .get_solid_id()] = cell;
                            }
                        }
                    }
                }
            }
        }
      return cells_cut_by_object;
    }
}


template typename DoFHandler<3>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(
  const DoFHandler<3> &dof_handler,
  const Point<3> &     point);
template typename DoFHandler<2>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_tree(
  const DoFHandler<2> &dof_handler,
  const Point<2> &     point);

template void
LetheGridTools::vertices_cell_mapping(
  const DoFHandler<2> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &vertices_cell_map);
template void
LetheGridTools::vertices_cell_mapping(
  const DoFHandler<3> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &vertices_cell_map);


template typename DoFHandler<2>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(
  const DoFHandler<2> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &                                                 vertices_cell_map,
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const Point<2> &                                    point);

template typename DoFHandler<3>::active_cell_iterator
LetheGridTools::find_cell_around_point_with_neighbors(
  const DoFHandler<3> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &                                                 vertices_cell_map,
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const Point<3> &                                    point);

template typename std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_cells_around_cell<2>(
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &                                                 vertices_cell_map,
  const typename DoFHandler<2>::active_cell_iterator &cell);

template typename std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_cells_around_cell<3>(
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &                                                 vertices_cell_map,
  const typename DoFHandler<3>::active_cell_iterator &cell);

template std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_cells_in_cells(
  const DoFHandler<2> &                               dof_handler,
  const typename DoFHandler<2>::active_cell_iterator &cell);

template std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_cells_in_cells(
  const DoFHandler<3> &                               dof_handler,
  const typename DoFHandler<3>::active_cell_iterator &cell);

template bool
LetheGridTools::cell_cut_by_flat<2>(
  const typename DoFHandler<2>::active_cell_iterator &       cell,
  const typename DoFHandler<2 - 1, 2>::active_cell_iterator &cell_flat);

template bool
LetheGridTools::cell_cut_by_flat<3>(
  const typename DoFHandler<3>::active_cell_iterator &       cell,
  const typename DoFHandler<3 - 1, 3>::active_cell_iterator &cell_flat);


template std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_cells_around_flat_cell(
  const DoFHandler<2> &                                      dof_handler,
  const typename DoFHandler<2 - 1, 2>::active_cell_iterator &cell,
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &vertices_cell_map);

template std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_cells_around_flat_cell(
  const DoFHandler<3> &                                      dof_handler,
  const typename DoFHandler<3 - 1, 3>::active_cell_iterator &cell,
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &vertices_cell_map);

template std::pair<std::pair<Point<2>, bool>, Tensor<1, 2>>
LetheGridTools::project_to_d_linear_object<2, 1>(
  const typename DoFHandler<1, 2>::active_cell_iterator &object,
  const Point<2> &                                       trial_point);

template std::pair<std::pair<Point<3>, bool>, Tensor<1, 3>>
LetheGridTools::project_to_d_linear_object<3, 2>(
  const typename DoFHandler<2, 3>::active_cell_iterator &object,
  const Point<3> &                                       trial_point);


template bool
LetheGridTools::cell_pierced_by_edge<3>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const TriaIterator<CellAccessor<1, 3>> &            cell_edge);

template bool
LetheGridTools::cell_pierced_by_edge<3>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  Point<3>                                            point_1,
  Point<3>                                            point_2);

template std::vector<typename DoFHandler<2>::active_cell_iterator>
LetheGridTools::find_boundary_cells_in_sphere(const DoFHandler<2> &dof_handler,
                                              const Point<2> &     center,
                                              const double         radius);
template std::vector<typename DoFHandler<3>::active_cell_iterator>
LetheGridTools::find_boundary_cells_in_sphere(const DoFHandler<3> &dof_handler,
                                              const Point<3> &     center,
                                              const double         radius);


template std::map<
  typename DoFHandler<3>::active_cell_iterator,
  std::map<unsigned int, typename DoFHandler<3, 3>::active_cell_iterator>>
LetheGridTools::find_cells_cut_by_object(
  const DoFHandler<3> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &                             vertices_cell_map,
  std::vector<SerialSolid<3, 3>> &list_of_objects);

template std::map<
  typename DoFHandler<2>::active_cell_iterator,
  std::map<unsigned int, typename DoFHandler<2, 2>::active_cell_iterator>>
LetheGridTools::find_cells_cut_by_object(
  const DoFHandler<2> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &                             vertices_cell_map,
  std::vector<SerialSolid<2, 2>> &list_of_objects);

template std::map<
  typename DoFHandler<3>::active_cell_iterator,
  std::map<unsigned int, typename DoFHandler<2, 3>::active_cell_iterator>>
LetheGridTools::find_cells_cut_by_object(
  const DoFHandler<3> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    &                             vertices_cell_map,
  std::vector<SerialSolid<2, 3>> &list_of_objects);

template std::map<
  typename DoFHandler<2>::active_cell_iterator,
  std::map<unsigned int, typename DoFHandler<1, 2>::active_cell_iterator>>
LetheGridTools::find_cells_cut_by_object(
  const DoFHandler<2> &dof_handler,
  std::map<unsigned int, std::set<typename DoFHandler<2>::active_cell_iterator>>
    &                             vertices_cell_map,
  std::vector<SerialSolid<1, 2>> &list_of_objects);


template <int dim>
std::tuple<std::vector<bool>, std::vector<Point<3>>, std::vector<Tensor<1, 3>>>
LetheGridTools::find_particle_triangle_projection(
  const std::vector<Point<dim>> &                      triangle,
  const std::vector<Particles::ParticleIterator<dim>> &particles,
  const unsigned int &                                 n_particles_in_base_cell)
{
  std::vector<bool>         pass_distance_check(n_particles_in_base_cell);
  std::vector<Point<3>>     projection_points(n_particles_in_base_cell);
  std::vector<Tensor<1, 3>> normal_vectors(n_particles_in_base_cell);

  auto &p_0 = triangle[0];
  auto &p_1 = triangle[1];
  auto &p_2 = triangle[2];

  const Tensor<1, dim> e_0 = p_1 - p_0;
  const Tensor<1, dim> e_1 = p_2 - p_0;

  const Tensor<1, dim> normal      = cross_product_3d(e_0, e_1);
  const double         norm_normal = normal.norm();
  Tensor<1, dim>       unit_normal = normal / norm_normal;
  Tensor<1, 3>         unit_normal_3d;
  Point<3>             pt_in_triangle_3d;

  const double a   = e_0.norm_square();
  const double b   = scalar_product(e_0, e_1);
  const double c   = e_1.norm_square();
  const double det = a * c - b * b;


  // Pre-allocation for speed
  Tensor<1, dim> vector_to_plane;
  Point<dim>     pt_in_triangle;

  unsigned int k = 0;
  for (auto &part : particles)
    {
      const double radius =
        part->get_properties()[DEM::PropertiesIndex::dp] * 0.5;
      Point<dim> particle_position = part->get_location();
      vector_to_plane              = p_0 - particle_position;

      // A bool variable for region 0
      bool region_zero = 0;

      // Check to see if the particle is located on the correct side (with
      // respect to the normal vector) of the triangle
      if (vector_to_plane * unit_normal > 0)
        {
          unit_normal = -1.0 * unit_normal;
        }

      double distance_squared = scalar_product(vector_to_plane, unit_normal);

      // If the particle is too far from the plane, set distance squared as an
      // arbitrary distance and continue
      if (distance_squared > (radius * radius))
        {
          pass_distance_check[k] = false;
          ++k;
          continue;
        }

      // Otherwise, do the full calculation taken from Eberly 2003
      const double d = scalar_product(e_0, vector_to_plane);
      const double e = scalar_product(e_1, vector_to_plane);

      // Calculate necessary values;
      double s = b * e - c * d;
      double t = b * d - a * e;

      if (s + t <= det)
        {
          if (s < 0)
            {
              if (t < 0)
                {
                  // Region 4
                  if (d < 0)
                    {
                      t = 0;
                      if (-d >= a)
                        s = 1;
                      else
                        s = -d / a;
                    }
                  else
                    {
                      s = 0;
                      if (e >= 0)
                        t = 0;
                      else if (-e >= c)
                        t = 1;
                      else
                        t = e / c;
                    }
                }
              else
                {
                  // Region 3
                  s = 0;
                  if (e >= 0)
                    t = 0;
                  else if (-e >= c)
                    t = 1;
                  else
                    t = -e / c;
                }
            }
          else if (t < 0)
            {
              // Region 5
              t = 0;
              if (d >= 0)
                s = 0;
              else if (-d >= a)
                s = 1;
              else
                s = -d / a;
            }
          else
            {
              // Region 0
              const double inv_det = 1. / det;
              s *= inv_det;
              t *= inv_det;

              // In region 0, normal vector is the face normal vector
              // Cast unit_normal to a tensor<1, 3>
              if constexpr (dim == 3)
                unit_normal_3d = unit_normal;

              if constexpr (dim == 2)
                unit_normal_3d = tensor_nd_to_3d(unit_normal);

              region_zero = 1;
            }
        }
      else
        {
          if (s < 0)
            {
              // Region 2
              const double tmp0 = b + d;
              const double tmp1 = c + e;
              if (tmp1 > tmp0)
                {
                  const double numer = tmp1 - tmp0;
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    s = 1;
                  else
                    s = numer / denom;

                  t = 1 - s;
                }
              else
                {
                  s = 0;
                  if (tmp1 <= 0)
                    t = 1;
                  else if (e >= 0)
                    t = 0;
                  else
                    t = -e / c;
                }
            }
          else if (t < 0)
            {
              // Region 6
              const double tmp0 = b + e;
              const double tmp1 = a + d;
              if (tmp1 > tmp0)
                {
                  const double numer = tmp1 - tmp0;
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    t = 1;
                  else
                    t = numer / denom;
                  s = 1 - t;
                }
              else
                {
                  t = 0;
                  if (tmp1 <= 0)
                    s = 1;
                  else if (d >= 0)
                    s = 0;
                  else
                    s = -d / a;
                }
            }
          else
            {
              // Region 1
              const double numer = (c + e) - (b + d);
              if (numer <= 0)
                s = 0;
              else
                {
                  const double denom = a - 2 * b + c;
                  if (numer >= denom)
                    s = 1;
                  else
                    s = numer / denom;
                }
              t = 1 - s;
            }
        }

      pt_in_triangle = p_0 + s * e_0 + t * e_1;

      if (!region_zero)
        {
          // normal vector
          const Tensor<1, dim> normal = particle_position - pt_in_triangle;

          if constexpr (dim == 3)
            unit_normal_3d = normal / normal.norm();

          if constexpr (dim == 2)
            unit_normal_3d = tensor_nd_to_3d(normal / normal.norm());
        }


      // Cast pt_in_triangle on Point<3>
      if constexpr (dim == 3)
        pt_in_triangle_3d = pt_in_triangle;

      if constexpr (dim == 2)
        pt_in_triangle_3d = point_nd_to_3d(pt_in_triangle);

      projection_points[k]   = pt_in_triangle_3d;
      pass_distance_check[k] = true;
      normal_vectors[k]      = unit_normal_3d;
      ++k;
    }
  return std::make_tuple(pass_distance_check,
                         projection_points,
                         normal_vectors);
}

template std::
  tuple<std::vector<bool>, std::vector<Point<3>>, std::vector<Tensor<1, 3>>>
  LetheGridTools::find_particle_triangle_projection(
    const std::vector<Point<2>> &                      triangle,
    const std::vector<Particles::ParticleIterator<2>> &particles,
    const unsigned int &n_particles_in_base_cell);
template std::
  tuple<std::vector<bool>, std::vector<Point<3>>, std::vector<Tensor<1, 3>>>
  LetheGridTools::find_particle_triangle_projection(
    const std::vector<Point<3>> &                      triangle,
    const std::vector<Particles::ParticleIterator<3>> &particles,
    const unsigned int &n_particles_in_base_cell);
