#include "dem/find_boundary_cells_information.h"

using namespace dealii;

// The constructor of this class is empty
template <int dim, int spacedim>
FindBoundaryCellsInformation<dim, spacedim>::FindBoundaryCellsInformation()
{}

// This function finds all the boundary cells and faces in the triangulation,
// for each cell the boundary faces are specified and the normal vector as well
// as a point on the boundary faces are obtained
template <int dim, int spacedim>
std::vector<boundary_cells_info_struct<dim>>
FindBoundaryCellsInformation<dim, spacedim>::find_boundary_cells_information(
  const Triangulation<dim, spacedim> &tr)
{
  // Initializing output_vector and a search_vector (containing boundary_id and
  // cell) to avoid replication of a boundary face. All the found boundary faces
  // will be searched in this vector before being added to the output_vector. If
  // they are not found in this search_vector they will be added to the
  // output_vector as well as the search_vector
  std::vector<boundary_cells_info_struct<dim>> output_vector;
  std::vector<std::pair<int, typename Triangulation<dim>::active_cell_iterator>>
    search_vector;

  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces
  const FE_Q<dim, spacedim> fe(1);
  QGauss<dim>               quadrature_formula(1);
  QGauss<dim - 1>           face_quadrature_formula(1);
  unsigned int              n_face_q_points = face_quadrature_formula.size();
  FEFaceValues<dim>         fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  // Iterating over the active cells in the trangulation
  for (typename Triangulation<dim>::active_cell_iterator cell_iterator =
         tr.begin_active();
       cell_iterator != tr.end();
       ++cell_iterator)
    {
      // Iterating over the faces of each cell
      for (unsigned int face_id = 0;
           face_id < GeometryInfo<dim>::faces_per_cell;
           ++face_id)
        {
          // Check to see if the face is located at boundary
          if (cell_iterator->face(face_id)->at_boundary())
            {
              fe_face_values.reinit(cell_iterator, face_id);

              for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                   ++f_q_point)
                {
                  // Finding the normal vector of the boundary face
                  Tensor<1, dim> normal_vector =
                    -1 * fe_face_values.normal_vector(f_q_point);

                  // Finding a point on the boundary face
                  Point<dim> quad_point = fe_face_values.quadrature_point(0);

                  // Boundary id
                  int boundary_id = cell_iterator->face(face_id)->boundary_id();

                  // Storing these information into the
                  // boundary_cells_info_struct
                  boundary_cells_info_struct<dim> boundary_information;
                  boundary_information.boundary_id      = boundary_id;
                  boundary_information.cell             = cell_iterator;
                  boundary_information.boundary_face_id = face_id;
                  boundary_information.normal_vector    = normal_vector;
                  boundary_information.point_on_face    = quad_point;

                  // Searching to see if boundary with these information were
                  // already found or not. If this is a new boundary face, it
                  // will be added to the output_vector as well as the
                  // search_vector to avoid repetition
                  auto information_search_element =
                    std::make_pair(boundary_id, cell_iterator);
                  auto search_iterator = std::find(search_vector.begin(),
                                                   search_vector.end(),
                                                   information_search_element);
                  if (search_iterator == search_vector.end())
                    {
                      output_vector.push_back(boundary_information);
                      search_vector.push_back(
                        std::make_pair(boundary_id, cell_iterator));
                    }
                }
            }
        }
    }
  return output_vector;
}

template class FindBoundaryCellsInformation<2, 2>;
template class FindBoundaryCellsInformation<3, 3>;
