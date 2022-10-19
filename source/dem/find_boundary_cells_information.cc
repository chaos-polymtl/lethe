#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/data_containers.h>
#include <dem/find_boundary_cells_information.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <boost/range/adaptor/map.hpp>


using namespace dealii;

// The constructor of this class is empty
template <int dim>
BoundaryCellsInformation<dim>::BoundaryCellsInformation()
{}

template <int dim>
void
BoundaryCellsInformation<dim>::build(
  const parallel::distributed::Triangulation<dim> & triangulation,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const std::vector<unsigned int> &                 outlet_boundaries,
  const bool &                                      check_diamond_cells,
  const bool &              expand_particle_wall_contact_search,
  const ConditionalOStream &pcout)
{
  boundary_cells_with_faces.clear();
  global_boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  local_cells_with_boundary_lines.clear();
  boundary_cells_information.clear();
  global_boundary_cells_information.clear();
  boundary_cells_with_points.clear();
  boundary_cells_for_floating_walls.clear();

  const double maximal_cell_diameter =
    GridTools::maximal_cell_diameter(triangulation);

  find_boundary_cells_information(triangulation, outlet_boundaries);

  // Finding boundary cells with lines and points
  find_particle_point_and_line_contact_cells(triangulation);

  // Finding cells adjacent to floating walls
  find_boundary_cells_for_floating_walls(triangulation,
                                         floating_wall_properties,
                                         maximal_cell_diameter);

  // Add cells with boundary lines to boundary cells container
  add_cells_with_boundary_lines_to_boundary_cells(triangulation,
                                                  outlet_boundaries,
                                                  check_diamond_cells,
                                                  pcout);

  if (expand_particle_wall_contact_search)
    {
      if (display_pw_contact_expansion_warning)
        {
          pcout
            << "Warning: expansion of particle-wall contact list is enabled. "
            << std::endl
            << "This feature should only be activated in geometries with concave boundaries. "
               "(For example, for particles flow inside a cylinder or sphere). In geometries with "
               "convex boundaries, this feature MUST NOT be activated."
            << std::endl;
          display_pw_contact_expansion_warning = false;
        }

      add_boundary_neighbors_of_boundary_cells(
        triangulation,
        outlet_boundaries,
        boundary_cells_information,
        global_boundary_cells_information);
    }
  else
    {
      if (display_pw_contact_expansion_warning)
        {
          pcout
            << "Warning: expansion of particle-wall contact list is disabled. "
            << std::endl
            << "This feature is useful in geometries with concave boundaries. "
            << std::endl;
          display_pw_contact_expansion_warning = false;
        }
    }
}

template <int dim>
void
BoundaryCellsInformation<dim>::build(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries,
  const bool &                                     check_diamond_cells,
  const ConditionalOStream &                       pcout)
{
  boundary_cells_with_faces.clear();
  global_boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  local_cells_with_boundary_lines.clear();
  boundary_cells_information.clear();
  global_boundary_cells_information.clear();
  boundary_cells_with_points.clear();
  boundary_cells_for_floating_walls.clear();

  find_boundary_cells_information(triangulation, outlet_boundaries);

  // Finding boundary cells with lines and points
  find_particle_point_and_line_contact_cells(triangulation);

  // Add cells with boundary lines to boundary cells container
  add_cells_with_boundary_lines_to_boundary_cells(triangulation,
                                                  outlet_boundaries,
                                                  check_diamond_cells,
                                                  pcout);
}

// This function finds all the boundary cells and faces in the triangulation,
// for each cell the boundary faces are specified and the normal vector as well
// as a point on the boundary faces are obtained
template <int dim>
void
BoundaryCellsInformation<dim>::find_boundary_cells_information(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries)
{
  // Initializing output_vector and a search_vector (containing boundary_id and
  // cell) to avoid replication of a boundary face. All the found boundary faces
  // will be searched in this vector before being added to the output_vector. If
  // they are not found in this search_vector they will be added to the
  // output_vector as well as the search_vector
  // std::map<int, boundary_cells_info_struct<dim>> output_vector;
  std::vector<std::pair<int, typename Triangulation<dim>::active_cell_iterator>>
    search_vector;

  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  // Iterating over the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // Iterating over the faces of each cell
      for (const auto &face : cell->face_iterators())
        {
          // We search to see if the boundary is defined as an outlet, periodic
          // or not. If it is not defined as one of those, we proceed.
          int  face_id = cell->face_iterator_to_index(face);
          bool is_outlet =
            std::find(outlet_boundaries.begin(),
                      outlet_boundaries.end(),
                      face->boundary_id()) != outlet_boundaries.end();
          bool is_periodic = cell->has_periodic_neighbor(face_id);

          if (!is_outlet && !is_periodic)
            {
              // Check if the face is located at boundary
              if (face->at_boundary())
                {
                  fe_face_values.reinit(cell, face_id);

                  // Finding the normal vector of the boundary face
                  Tensor<1, dim> normal_vector =
                    -fe_face_values.normal_vector(0);

                  // Finding a point on the boundary face
                  Point<dim> quad_point = fe_face_values.quadrature_point(0);

                  // Get global index of the face
                  unsigned int global_face_id = cell->face_index(face_id);

                  // Storing these information into the
                  // boundary_cells_info_struct
                  boundary_cells_info_struct<dim> boundary_information;
                  boundary_information.cell           = cell;
                  boundary_information.boundary_id    = face->boundary_id();
                  boundary_information.global_face_id = global_face_id;
                  boundary_information.normal_vector  = normal_vector;
                  boundary_information.point_on_face  = quad_point;

                  // Searching to see if boundary with this information was
                  // already found or not. If this is a new
                  // boundary face, it will be added to the
                  // output_vector as well as the search_vector to avoid
                  // repetition
                  auto information_search_element =
                    std::make_pair(face_id, cell);
                  auto search_iterator = std::find(search_vector.begin(),
                                                   search_vector.end(),
                                                   information_search_element);
                  if (search_iterator == search_vector.end())
                    {
                      if (cell->is_locally_owned())
                        {
                          boundary_cells_information.insert(
                            {global_face_id, boundary_information});
                          boundary_cells_with_faces.push_back(cell);
                          search_vector.push_back(information_search_element);
                        }
                      global_boundary_cells_information.insert(
                        {global_face_id, boundary_information});
                    }
                }
            }
        }
    }
}


// This function is used to update the normal vector and the location of the
// stored point of the boundary faces. It is used when the grid is moving.
// Updated points and normal vectors are then used to update the particle-wall
// contact list.
template <int dim>
void
BoundaryCellsInformation<dim>::update_boundary_info_after_grid_motion(
  typename DEM::dem_data_structures<dim>::boundary_points_and_normal_vectors
    &updated_boundary_points_and_normal_vectors)
{
  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);
  for (auto &cell : boundary_cells_with_faces)
    {
      if (cell->is_locally_owned())
        {
          // Iterating over the faces of each cell
          for (const auto &face : cell->face_iterators())
            {
              // Check to see if the face is located at boundary
              if (face->at_boundary())
                {
                  int face_id = cell->face_iterator_to_index(face);
                  fe_face_values.reinit(cell, face_id);

                  // Finding the normal vector of the boundary face
                  Tensor<1, 3> normal_vector;

                  if constexpr (dim == 3)
                    normal_vector = -fe_face_values.normal_vector(0);

                  if constexpr (dim == 2)
                    normal_vector =
                      tensor_nd_to_3d(-fe_face_values.normal_vector(0));

                  // Finding a point on the boundary face
                  Point<3> quad_point;

                  if constexpr (dim == 3)
                    quad_point = fe_face_values.quadrature_point(0);

                  if constexpr (dim == 2)
                    quad_point =
                      point_nd_to_3d(fe_face_values.quadrature_point(0));

                  updated_boundary_points_and_normal_vectors[cell->face_index(
                    face_id)] = std::make_pair(normal_vector, quad_point);
                }
            }
        }
    }
}

// This function returns all the cells at which particle-line and particle-point
// contact search should be performed
template <int dim>
void
BoundaryCellsInformation<dim>::find_particle_point_and_line_contact_cells(
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  // Here, we first loop through all the non-boundary faces of the local cells
  // and find (add them to all cells_with_boundary_lines container) all boundary
  // lines of these faces. Then the boundary lines which are the borders of two
  // boundary faces should be removed from this container. Finally the remained
  // elements of this container are added to boundary_cells_with_lines

  // Boundary lines only exist in three-dimensional cases
  if (dim == 3)
    {
      std::unordered_map<
        std::string,
        std::unordered_map<unsigned int, std::pair<Point<dim>, Point<dim>>>>
        all_cells_with_boundary_lines;

      // Iterating over the active cells in the triangulation
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              // Looping through all the faces of the cell
              for (const auto &face : cell->face_iterators())
                {
                  // Check to see if the face is not located at boundary
                  if (!face->at_boundary())
                    {
                      // Adding all the boundary lines of these faces into the
                      // all_cells_with_boundary_lines container
                      for (unsigned int l = 0; l < face->n_lines(); ++l)
                        {
                          if (face->line(l)->at_boundary())
                            {
                              all_cells_with_boundary_lines[cell->id()
                                                              .to_string()]
                                .insert(
                                  {face->line_index(l),
                                   std::make_pair(face->line(l)->vertex(0),
                                                  face->line(l)->vertex(1))});
                            }
                        }
                    }
                }
            }
        }

      // Now iterating again over the cells, faces and lines and if the boundary
      // line of the boundary face exists in the container, it is removed from
      // the container
      if (!all_cells_with_boundary_lines.empty())
        {
          for (const auto &cell : triangulation.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  for (const auto &face : cell->face_iterators())
                    {
                      // Check to see if the face is located at boundary
                      if (face->at_boundary())
                        {
                          for (unsigned int l = 0; l < face->n_lines(); ++l)
                            {
                              all_cells_with_boundary_lines[cell->id()
                                                              .to_string()]
                                .erase(face->line_index(l));
                            }
                        }
                    }
                }
            }
        }

      // Finally, adding the remained elements of the
      // all_cells_with_boundary_lines container to boundary_cells_with_lines
      if (!all_cells_with_boundary_lines.empty())
        {
          for (const auto &cell : triangulation.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  std::string cell_id_string = cell->id().to_string();
                  std::unordered_map<unsigned int,
                                     std::pair<Point<dim>, Point<dim>>>
                    &cell_boundary_lines =
                      all_cells_with_boundary_lines[cell_id_string];

                  if (!cell_boundary_lines.empty())
                    {
                      for (auto &map_iterator : cell_boundary_lines)
                        {
                          boundary_cells_with_lines.insert(
                            {cell_id_string,
                             std::make_tuple(cell,
                                             map_iterator.second.first,
                                             map_iterator.second.second)});

                          // Add the cell to local_cells_with_boundary_lines to
                          // be used in
                          // add_cells_with_boundary_lines_to_boundary_cells
                          local_cells_with_boundary_lines.push_back(cell);
                        }
                    }
                }
            }
        }
    }

  // Finding boundary points. We need this container because at_boundary()
  // function is not defined for vertex id (unsigned int) nor vertex position
  // (Point<dim>). For line and face objects, at_boundary() function is
  // available to help

  // First getting a set of all boundary vertices
  std::set<unsigned int> boundary_vertices;
  for (const auto &face : triangulation.active_face_iterators())
    {
      if (face->at_boundary())
        {
          for (unsigned int v = 0; v < face->n_vertices(); ++v)
            {
              boundary_vertices.insert(face->vertex_index(v));
            }
        }
    }

  // If a cell does not have any boundary line, but has boundary vertex, it is
  // added to boundary_cells_with_points container Iterating over the active
  // cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() && !cell->has_boundary_lines())
        {
          for (const auto &face : cell->face_iterators())
            {
              if (!face->at_boundary())
                {
                  for (unsigned int v = 0; v < face->n_vertices(); ++v)
                    {
                      if (boundary_vertices.count(face->vertex_index(v)) > 0)
                        {
                          boundary_cells_with_points.insert(
                            {cell->id().to_string(),
                             std::make_pair(cell, face->vertex(v))});
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
BoundaryCellsInformation<dim>::find_global_boundary_cells_with_faces(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries)
{
  // Iterating over the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // Iterating over the faces of each cell
      for (const auto &face : cell->face_iterators())
        {
          // We search to see if the boundary is defined as an outlet or
          // not. If it is not defined as an outlet we proceed.
          if (std::find(outlet_boundaries.begin(),
                        outlet_boundaries.end(),
                        face->boundary_id()) == outlet_boundaries.end())
            {
              // Check to see if the face is located at boundary
              if (face->at_boundary())
                {
                  global_boundary_cells_with_faces.push_back(cell);
                }
            }
        }
    }
}

template <int dim>
void
BoundaryCellsInformation<dim>::add_cells_with_boundary_lines_to_boundary_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries,
  const bool &                                     check_diamond_cells,
  const ConditionalOStream &                       pcout)
{
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    boundary_neighbor_cells;

  // First we have to store all the global (including local and non-local) cells
  // with boundary faces in a container. Note that the boundary_cells_with_faces
  // container only stores the local boundary cells
  find_global_boundary_cells_with_faces(triangulation, outlet_boundaries);

  if (!local_cells_with_boundary_lines.empty())
    {
      for (auto &cell_with_boundary_line : local_cells_with_boundary_lines)
        {
          // Each cell with boundary line must have two neighbor cells with
          // boundary faces. We have to find these two neighbors first.
          // Finding neighbors of the cell
          std::vector<typename Triangulation<dim>::active_cell_iterator>
            active_neighbors;
          GridTools::get_active_neighbors<Triangulation<dim>>(
            cell_with_boundary_line, active_neighbors);

          // Looping through the neighbors and storing the ones with boundary
          // faces
          boundary_neighbor_cells.clear();
          for (auto &neighbor : active_neighbors)
            {
              if (std::count(global_boundary_cells_with_faces.begin(),
                             global_boundary_cells_with_faces.end(),
                             neighbor))
                {
                  boundary_neighbor_cells.push_back(neighbor);
                }
            }

          // Find the information of every combination of two cells in
          // global_boundary_cells_normals with counters
          for (unsigned int counter_one = 0;
               counter_one < boundary_neighbor_cells.size();
               ++counter_one)
            {
              for (unsigned int counter_two = 0;
                   counter_two < boundary_neighbor_cells.size();
                   ++counter_two)
                {
                  if (counter_one > counter_two)
                    {
                      // Get the pair of cells
                      typename Triangulation<dim>::active_cell_iterator
                        cell_one = boundary_neighbor_cells[counter_one];
                      typename Triangulation<dim>::active_cell_iterator
                        cell_two = boundary_neighbor_cells[counter_two];

                      // Loop on all faces of the 2 cells
                      for (const auto &face_one : cell_one->face_iterators())
                        {
                          for (const auto &face_two :
                               cell_two->face_iterators())
                            {
                              // Check if faces are located at boundaries
                              if (face_one->at_boundary() &&
                                  face_two->at_boundary())
                                {
                                  unsigned int face_id_one, face_id_two;
                                  face_id_one =
                                    cell_one->face_iterator_to_index(face_one);
                                  face_id_two =
                                    cell_two->face_iterator_to_index(face_two);

                                  // Check if faces are on a wall boundary
                                  // (not outlet nor periodic)
                                  bool face_one_is_wall, face_two_is_wall;
                                  face_one_is_wall =
                                    (std::find(outlet_boundaries.begin(),
                                               outlet_boundaries.end(),
                                               face_one->boundary_id()) ==
                                       outlet_boundaries.end() &&
                                     !cell_one->has_periodic_neighbor(
                                       face_id_one));
                                  face_two_is_wall =
                                    (std::find(outlet_boundaries.begin(),
                                               outlet_boundaries.end(),
                                               face_two->boundary_id()) ==
                                       outlet_boundaries.end() &&
                                     !cell_two->has_periodic_neighbor(
                                       face_id_two));

                                  if (face_one_is_wall && face_two_is_wall)
                                    {
                                      // Get the normal vector of the first cell
                                      Tensor<1, dim> normal_vector_one,
                                        normal_vector_two;
                                      normal_vector_one =
                                        global_boundary_cells_information
                                          .at(cell_one->face_index(face_id_one))
                                          .normal_vector;

                                      normal_vector_two =
                                        global_boundary_cells_information
                                          .at(cell_two->face_index(face_id_two))
                                          .normal_vector;

                                      // Check to see if the dot product of the
                                      // two normal vectors is larger than a
                                      // specified threshold (cos(45) = 0.707)
                                      if (normal_vector_one *
                                            normal_vector_two >
                                          0.707)
                                        {
                                          if (display_diamond_cells_warning)
                                            {
                                              pcout
                                                << std::endl
                                                << "Warning: There are diamond-shaped cells in the input triangulation. It is strongly recommended to use a different triangulation without such cells. It should be mentioned that these cells are not detected if you have grid motion"
                                                << std::endl;
                                              display_diamond_cells_warning =
                                                false;
                                            }

                                          if (check_diamond_cells &&
                                              cell_with_boundary_line
                                                ->is_locally_owned())
                                            {
                                              // If this condition is true, add
                                              // this cell to
                                              // boundary_cells_information.
                                              // Since the key of
                                              // boundary_cells_information is
                                              // the boundary face id, and this
                                              // cell does not have a boundary
                                              // face, we add it to the
                                              // boundary_cells_information with
                                              // a negative key
                                              int imaginary_face_id =
                                                -1 * cell_two->face_index(
                                                       face_id_two);

                                              // Update cell to the cell with
                                              // boundary line
                                              boundary_cells_info_struct<dim>
                                                imaginary_information =
                                                  global_boundary_cells_information
                                                    .at(cell_two->face_index(
                                                      face_id_two));
                                              imaginary_information.cell =
                                                cell_with_boundary_line;
                                              imaginary_information
                                                .global_face_id =
                                                imaginary_face_id;

                                              boundary_cells_information.insert(
                                                {imaginary_face_id,
                                                 imaginary_information});
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// This function finds the triangulation cells adjacent to the floating walls
template <int dim>
void
BoundaryCellsInformation<dim>::find_boundary_cells_for_floating_walls(
  const parallel::distributed::Triangulation<dim> & triangulation,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      maximum_cell_diameter)
{
  // Reading floating wall properties
  std::vector<Point<dim>> point_on_wall =
    floating_wall_properties.points_on_walls;
  std::vector<Tensor<1, dim>> wall_normal_vector =
    floating_wall_properties.floating_walls_normal_vectors;
  unsigned int floating_wall_number =
    floating_wall_properties.floating_walls_number;

  // Looping through floating walls
  for (unsigned int floating_wall_id = 0;
       floating_wall_id < floating_wall_number;
       ++floating_wall_id)
    {
      // Looping through cells
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          // If the cell is owned by the processor
          if (cell->is_locally_owned())
            {
              // Looping over vertices of each cell
              for (unsigned int vertex = 0; vertex < cell->n_vertices();
                   ++vertex)
                {
                  // Finding vertex-floating wall distance
                  Tensor<1, dim> connecting_vector =
                    cell->vertex(vertex) - point_on_wall[floating_wall_id];
                  double vertex_wall_distance =
                    connecting_vector * wall_normal_vector[floating_wall_id];

                  // If the distance is less than the largest cell size,
                  // it should be added to the
                  // boundary_cells_for_floating_walls
                  if (abs(vertex_wall_distance) < maximum_cell_diameter)
                    {
                      boundary_cells_for_floating_walls[floating_wall_id]
                        .insert(cell);
                    }
                }
            }
        }
    }
}

template <int dim>
void
BoundaryCellsInformation<dim>::add_boundary_neighbors_of_boundary_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries,
  std::map<int, boundary_cells_info_struct<dim>> & boundary_cells_information,
  const std::map<int, boundary_cells_info_struct<dim>>
    &global_boundary_cells_information)
{
  // Create a vector of a set of adjacent cells of all the vertices
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    v_to_c = GridTools::vertex_to_cell_map(triangulation);

  for (auto &&boundary_cells_info :
       boundary_cells_information | boost::adaptors::map_values)
    {
      // The boundary cell is local by definition
      // This loop will check for all the cell neighbors located at vertices of
      // a boundary cell and will add the proper information in the boundary
      // cell list if faces are not at outlet nor periodic boundary.

      // Iterate over the vertices of each boundary cell
      for (auto vertex_id : GeometryInfo<dim>::vertex_indices())
        {
          // Iterate over the neighbors of each boundary cell having vertex_id
          for (const auto &neighbor :
               v_to_c[boundary_cells_info.cell->vertex_index(vertex_id)])
            {
              // Iterate over the faces of the neighbor cell
              for (const auto &face : neighbor->face_iterators())
                {
                  // Check if the face of the neighbor is located at boundary
                  if (face->at_boundary())
                    {
                      // Check if face is on a wall boundary (not outlet nor
                      // periodic)
                      unsigned int face_id =
                        neighbor->face_iterator_to_index(face);
                      bool face_is_wall =
                        (std::find(outlet_boundaries.begin(),
                                   outlet_boundaries.end(),
                                   face->boundary_id()) ==
                           outlet_boundaries.end() &&
                         !neighbor->has_periodic_neighbor(face_id));

                      if (face_is_wall)
                        {
                          // Get the boundary neighbor cell info
                          boundary_cells_info_struct<dim>
                            boundary_neighbor_information =
                              global_boundary_cells_information.at(
                                neighbor->face_index(face_id));

                          // Add the main boundary cell with the information
                          // (point and normal vector) of the neighbor boundary
                          // cell to the boundary_cells_information container.
                          // Note that since we may already have an element with
                          // the key of face_id (key of the
                          // boundary_cells_information map) in the
                          // boundary_cells_information, we add the new element
                          // with a unique key to create a unique id in the map.
                          // This unique key is generated using Cantor pairing
                          // function: unique_key = 0.5 * (a + b) * (a + b + 1)
                          // + b where a and b are global boundary face ids of
                          // the main boundary cell and the neighbor boundary
                          // cell.
                          int imaginary_face_id =
                            -0.5 *
                              (boundary_cells_info.global_face_id +
                               boundary_neighbor_information.global_face_id) *
                              (boundary_cells_info.global_face_id +
                               boundary_neighbor_information.global_face_id +
                               1) +
                            boundary_neighbor_information.global_face_id;

                          // Create a cell info object which is a copy of all
                          // the boundary neighbor information applied to the
                          // boundary cell & store in map with imaginary key.
                          boundary_cells_info_struct<dim> boundary_information =
                            boundary_neighbor_information;
                          boundary_information.cell = boundary_cells_info.cell;

                          boundary_cells_information.insert(
                            {imaginary_face_id, boundary_information});
                        }
                    }
                }
            }
        }
    }
}

template class BoundaryCellsInformation<2>;
template class BoundaryCellsInformation<3>;
