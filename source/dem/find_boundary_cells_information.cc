#include <dem/find_boundary_cells_information.h>

using namespace dealii;

// The constructor of this class is empty
template <int dim>
FindBoundaryCellsInformation<dim>::FindBoundaryCellsInformation()
{}

// This function finds all the boundary cells and faces in the triangulation,
// for each cell the boundary faces are specified and the normal vector as well
// as a point on the boundary faces are obtained
template <int dim>
std::map<int, boundary_cells_info_struct<dim>>
FindBoundaryCellsInformation<dim>::find_boundary_cells_information(
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    &                                              boundary_cells_with_faces,
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  // Initializing output_vector and a search_vector (containing boundary_id and
  // cell) to avoid replication of a boundary face. All the found boundary faces
  // will be searched in this vector before being added to the output_vector. If
  // they are not found in this search_vector they will be added to the
  // output_vector as well as the search_vector
  std::map<int, boundary_cells_info_struct<dim>> output_vector;
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
                        -1 * fe_face_values.normal_vector(f_q_point);

                      // Finding a point on the boundary face
                      Point<dim> quad_point =
                        fe_face_values.quadrature_point(0);

                      // Storing these information into the
                      // boundary_cells_info_struct
                      boundary_cells_info_struct<dim> boundary_information;
                      boundary_information.cell = cell;
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
                          output_vector.insert(
                            {cell->face_index(face_id), boundary_information});
                          boundary_cells_with_faces.push_back(cell);
                          search_vector.push_back(information_search_element);
                        }
                    }
                }
            }
        }
    }
  return output_vector;
}

// This function returns all the cells at which particle-line and particle-point
// contact search should be performed
template <int dim>
void
FindBoundaryCellsInformation<dim>::find_particle_point_and_line_contact_cells(
  const std::vector<typename Triangulation<dim>::active_cell_iterator>
    &                                              boundary_cells_with_faces,
  const parallel::distributed::Triangulation<dim> &triangulation,
  std::vector<std::tuple<typename Triangulation<dim>::active_cell_iterator,
                         Point<dim>,
                         Point<dim>>> &            boundary_cells_with_lines,
  std::vector<std::pair<typename Triangulation<dim>::active_cell_iterator,
                        Point<dim>>> &             boundary_cells_with_points)
{
  // Looing over all the faces to find boundary faces and then looping over the
  // vertices of these boundary faces to find all the vertices located on
  // boundaries
  for (const auto &face : triangulation.active_face_iterators())
    {
      if (face->at_boundary())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
               ++v)
            {
              boundary_vertices.insert(
                {face->vertex_index(v), face->vertex(v)});
            }
        }
    }

  // Looping over boundary vertices, finding all the adjacent cells to these
  // vertices and if these cells (which atleast own one boundary vertex) do not
  // exist in the boundary_cells_with_faces, they will be stored in
  // boundary_cells_with_lines_or_points
  for (auto iterator = boundary_vertices.begin();
       iterator != boundary_vertices.end();
       ++iterator)
    {
      auto vertex_index = iterator->first;
      auto candidate_cells =
        GridTools::find_cells_adjacent_to_vertex(triangulation, vertex_index);
      for (unsigned int counter = 0; counter != candidate_cells.size();
           ++counter)
        {
          auto search_iterator = std::find(boundary_cells_with_faces.begin(),
                                           boundary_cells_with_faces.end(),
                                           candidate_cells[counter]);
          if (search_iterator == boundary_cells_with_faces.end())
            {
              boundary_cells_with_lines_or_points.push_back(
                candidate_cells[counter]);
            }
        }
    }

  // Looping over boundary_cells_with_lines_or_points and counting the number of
  // boundary vertices for each cell. If the cell have one boundary vertex it
  // will be stored in boundary_cells_with_points, and if it has two boundary
  // vertices in boundary_cells_with_lines
  // The location of these boundary vertices are also stored to be used for
  // contact detection (fine search)
  for (unsigned int counter = 0;
       counter != boundary_cells_with_lines_or_points.size();
       ++counter)
    {
      number_of_boundary_vertices = 0;
      boundary_points.clear();
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
          boundary_cells_with_points.push_back(
            std::make_pair(cell, *boundary_points.begin()));
        }
      else if (number_of_boundary_vertices == 2)
        {
          boundary_cells_with_lines.push_back(std::make_tuple(
            cell, *boundary_points.begin(), *boundary_points.end()));
        }
    }
}

template class FindBoundaryCellsInformation<2>;
template class FindBoundaryCellsInformation<3>;
