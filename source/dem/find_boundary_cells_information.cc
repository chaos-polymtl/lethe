#include <dem/find_boundary_cells_information.h>

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
  const std::vector<unsigned int> &                 outlet_boundaries)
{
  boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  boundary_cells_information.clear();
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
}

template <int dim>
void
BoundaryCellsInformation<dim>::build(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<unsigned int> &                outlet_boundaries)
{
  boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  boundary_cells_information.clear();
  boundary_cells_with_points.clear();
  boundary_cells_for_floating_walls.clear();

  find_boundary_cells_information(triangulation, outlet_boundaries);

  // Finding boundary cells with lines and points
  find_particle_point_and_line_contact_cells(triangulation);
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
  unsigned int      n_face_q_points = face_quadrature_formula.size();
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Iterating over the faces of each cell
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              // We search to see if the boundary is defined as an outlet or
              // not. If it is not defined as an outlet we proceed.
              if (std::find(outlet_boundaries.begin(),
                            outlet_boundaries.end(),
                            cell->face(face_id)->boundary_id()) ==
                  outlet_boundaries.end())
                {
                  // Check to see if the face is located at boundary
                  if (cell->face(face_id)->at_boundary())
                    {
                      fe_face_values.reinit(cell, face_id);

                      for (unsigned int f_q_point = 0;
                           f_q_point < n_face_q_points;
                           ++f_q_point)
                        {
                          // Finding the normal vector of the boundary face
                          Tensor<1, dim> normal_vector =
                            -fe_face_values.normal_vector(f_q_point);

                          // Finding a point on the boundary face
                          Point<dim> quad_point =
                            fe_face_values.quadrature_point(0);

                          // Storing these information into the
                          // boundary_cells_info_struct
                          boundary_cells_info_struct<dim> boundary_information;
                          boundary_information.cell        = cell;
                          boundary_information.boundary_id = cell->face(face_id)->boundary_id();
                          boundary_information.global_face_id =
                            cell->face_index(face_id);
                          boundary_information.normal_vector = normal_vector;
                          boundary_information.point_on_face = quad_point;

                          // Searching to see if boundary with these information
                          // were already found or not. If this is a new
                          // boundary face, it will be added to the
                          // output_vector as well as the search_vector to avoid
                          // repetition
                          auto information_search_element =
                            std::make_pair(face_id, cell);
                          auto search_iterator =
                            std::find(search_vector.begin(),
                                      search_vector.end(),
                                      information_search_element);
                          if (search_iterator == search_vector.end())
                            {
                              boundary_cells_information.insert(
                                {cell->face_index(face_id),
                                 boundary_information});
                              boundary_cells_with_faces.push_back(cell);
                              search_vector.push_back(
                                information_search_element);
                            }
                        }
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
  std::map<unsigned int, std::pair<Tensor<1, dim>, Point<dim>>>
    &updated_boundary_points_and_normal_vectors)
{
  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  unsigned int      n_face_q_points = face_quadrature_formula.size();
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  for (auto &cell : boundary_cells_with_faces)
    {
      // Iterating over the faces of each cell
      for (int face_id = 0; face_id < int(GeometryInfo<dim>::faces_per_cell);
           ++face_id)
        {
          // Check to see if the face is located at boundary
          if (cell->face(face_id)->at_boundary())
            {
              fe_face_values.reinit(cell, face_id);

              for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                   ++f_q_point)
                {
                  // Finding the normal vector of the boundary face
                  Tensor<1, dim> normal_vector =
                    -fe_face_values.normal_vector(f_q_point);

                  // Finding a point on the boundary face
                  Point<dim> quad_point = fe_face_values.quadrature_point(0);
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

      // Iterating over the active cells in the trangulation
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              // Looping through all the faces of the cell
              for (int face_id = 0;
                   face_id < int(GeometryInfo<dim>::faces_per_cell);
                   ++face_id)
                {
                  // Check to see if the face is not located at boundary
                  if (!(cell->face(face_id)->at_boundary()))
                    {
                      // Adding all the boundary lines of these faces into the
                      // all_cells_with_boundary_lines container
                      for (unsigned int l = 0;
                           l < GeometryInfo<dim>::lines_per_face;
                           ++l)
                        {
                          if (cell->face(face_id)->line(l)->at_boundary())
                            {
                              all_cells_with_boundary_lines[cell->id()
                                                              .to_string()]
                                .insert(
                                  {cell->face(face_id)->line_index(l),
                                   std::make_pair(
                                     cell->face(face_id)->line(l)->vertex(0),
                                     cell->face(face_id)->line(l)->vertex(1))});
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
                  for (int face_id = 0;
                       face_id < int(GeometryInfo<dim>::faces_per_cell);
                       ++face_id)
                    {
                      // Check to see if the face is located at boundary
                      if ((cell->face(face_id)->at_boundary()))
                        {
                          for (unsigned int l = 0;
                               l < GeometryInfo<dim>::lines_per_face;
                               ++l)
                            {
                              all_cells_with_boundary_lines[cell->id()
                                                              .to_string()]
                                .erase(cell->face(face_id)->line_index(l));
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
                  auto &      cell_boundary_lines =
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
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
               ++v)
            {
              boundary_vertices.insert(face->vertex_index(v));
            }
        }
    }

  // If a cell does not have any boundary line, but has boundary vertex, it is
  // added to boundary_cells_with_points container Iterating over the active
  // cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() && !(cell->has_boundary_lines()))
        {
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              if (!(cell->face(face_id)->at_boundary()))
                {
                  for (unsigned int v = 0;
                       v < GeometryInfo<dim>::vertices_per_face;
                       ++v)
                    {
                      if (boundary_vertices.count(
                            cell->face(face_id)->vertex_index(v)) > 0)
                        {
                          boundary_cells_with_points.insert(
                            {cell->id().to_string(),
                             std::make_pair(cell,
                                            cell->face(face_id)->vertex(v))});
                        }
                    }
                }
            }
        }
    }
}

// This function finds the triangulation cells adjacent to the floating
// walls
template <int dim>
void
BoundaryCellsInformation<dim>::find_boundary_cells_for_floating_walls(
  const parallel::distributed::Triangulation<dim> & triangulation,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double &                                    maximum_cell_diameter)
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
              for (unsigned int vertex = 0;
                   vertex < GeometryInfo<dim>::vertices_per_cell;
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

template class BoundaryCellsInformation<2>;
template class BoundaryCellsInformation<3>;
