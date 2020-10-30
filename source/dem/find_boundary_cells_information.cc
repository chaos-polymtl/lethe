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
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties)
{
  boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  boundary_cells_information.clear();
  boundary_cells_with_points.clear();
  boundary_cells_for_floating_walls.clear();

  const double maximal_cell_diameter =
    GridTools::maximal_cell_diameter(triangulation);

  find_boundary_cells_information(triangulation);

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
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  boundary_cells_with_faces.clear();
  boundary_cells_with_lines.clear();
  boundary_cells_information.clear();
  boundary_cells_with_points.clear();
  boundary_cells_for_floating_walls.clear();

  find_boundary_cells_information(triangulation);

  // Finding boundary cells with lines and points
  find_particle_point_and_line_contact_cells(triangulation);
}

// This function finds all the boundary cells and faces in the triangulation,
// for each cell the boundary faces are specified and the normal vector as well
// as a point on the boundary faces are obtained
template <int dim>
void
BoundaryCellsInformation<dim>::find_boundary_cells_information(
  const parallel::distributed::Triangulation<dim> &triangulation)
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
                      Point<dim> quad_point =
                        fe_face_values.quadrature_point(0);

                      // Storing these information into the
                      // boundary_cells_info_struct
                      boundary_cells_info_struct<dim> boundary_information;
                      boundary_information.cell        = cell;
                      boundary_information.boundary_id = face_id;
                      boundary_information.boundary_face_id =
                        cell->face_index(face_id);
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_face = quad_point;

                      // Searching to see if boundary with these information
                      // were already found or not. If this is a new boundary
                      // face, it will be added to the output_vector as well as
                      // the search_vector to avoid repetition
                      auto information_search_element =
                        std::make_pair(face_id, cell);
                      auto search_iterator =
                        std::find(search_vector.begin(),
                                  search_vector.end(),
                                  information_search_element);
                      if (search_iterator == search_vector.end())
                        {
                          boundary_cells_information.insert(
                            {cell->face_index(face_id), boundary_information});
                          boundary_cells_with_faces.push_back(cell);
                          search_vector.push_back(information_search_element);
                        }
                    }
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
  // This vector stores both the cells with boundary lines and cells with
  // boundary points
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    boundary_cells_with_lines_or_points;

  // This unordered map stores the vertex index and position of boundary
  // vertices
  std::unordered_map<int, Point<dim>> boundary_vertices;

  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          // Iterating over the faces of each cell
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              // Check to see if the face is located at boundary
              if (cell->face(face_id)->at_boundary())
                {
                  for (unsigned int v = 0;
                       v < GeometryInfo<dim>::vertices_per_face;
                       ++v)
                    {
                      boundary_vertices.insert(
                        {cell->face(face_id)->vertex_index(v),
                         cell->face(face_id)->vertex(v)});
                    }
                }
            }
        }
    }

  // Looping over boundary vertices, finding all the adjacent cells to these
  // vertices and if these cells (which atleast own one boundary vertex) do not
  // exist in the boundary_cells_with_faces, they will be stored in
  // boundary_cells_with_lines_or_points

  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  for (auto iterator = boundary_vertices.begin();
       iterator != boundary_vertices.end();
       ++iterator)
    {
      auto vertex_index = iterator->first;

      for (const auto &neighbor : v_to_c[vertex_index])
        {
          if (neighbor->is_locally_owned() && !neighbor->at_boundary())
            {
              boundary_cells_with_lines_or_points.push_back(neighbor);
            }
        }
    }

  // Looping over boundary_cells_with_lines_or_points and counting the
  // number of boundary vertices for each cell. If the cell have one
  // boundary vertex it will be stored in boundary_cells_with_points,
  // and if it has two boundary vertices in boundary_cells_with_lines
  // The location of these boundary vertices are also stored to be used
  // for contact detection (fine search)
  for (unsigned int counter = 0;
       counter != boundary_cells_with_lines_or_points.size();
       ++counter)
    {
      unsigned int number_of_boundary_vertices = 0;

      // This vector stores the location of vertices on boundaries for
      // each cell. The size of this vector can be 1 or 2, since cells
      // with points have one boundary vertex and cells with lines have
      // two boundary vertices
      std::vector<Point<dim>> boundary_points;
      auto cell = boundary_cells_with_lines_or_points[counter];
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          if (boundary_vertices.count(cell->vertex_index(v)) > 0)
            {
              boundary_points.push_back(cell->vertex(v));
              ++number_of_boundary_vertices;
            }
        }
      if (number_of_boundary_vertices == 1)
        {
          // It means the cells has a boundary point
          boundary_cells_with_points.push_back(
            std::make_pair(cell, *boundary_points.begin()));
        }
      else if (number_of_boundary_vertices == 2)
        {
          // It means the cells has a boundary line
          boundary_cells_with_lines.push_back(std::make_tuple(
            cell, *boundary_points.begin(), *boundary_points.end()));
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
