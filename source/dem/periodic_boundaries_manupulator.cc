#include <core/tensors_and_points_dimension_manipulation.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <dem/periodic_boundaries_manupulator.h>

using namespace dealii;

template <int dim>
PeriodicBoundariesManipulator<dim>::PeriodicBoundariesManipulator()
{}



template <int dim>
void
PeriodicBoundariesManipulator<dim>::map_periodic_cells(
  const parallel::distributed::Triangulation<dim> & triangulation,
  const std::vector<unsigned int> &                 outlet_boundaries)
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

  // In case periodic faces
  FEFaceValues<dim> fe_face_values_periodic(fe,
                                            face_quadrature_formula,
                                            update_values |
                                              update_quadrature_points |
                                              update_normal_vectors);
  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // Iterating over the faces of each cell
      for (int face_id = 0; face_id < int(GeometryInfo<dim>::faces_per_cell);
           ++face_id)
        {
          // We search to see if the boundary is defined as an outlet or
          // not. If it is not defined as an outlet we proceed.
          if (cell->face(face_id)->at_boundary())
            {
              bool boundary_is_outlet =
                std::find(outlet_boundaries.begin(),
                          outlet_boundaries.end(),
                          cell->face(face_id)->boundary_id()) !=
                outlet_boundaries.end();

              if (cell->has_periodic_neighbor(face_id))
                {
                  typename Triangulation<dim>::active_cell_iterator
                    periodic_cell = cell->periodic_neighbor(face_id);

                  int periodic_face = cell->periodic_neighbor_face_no(face_id);

                  // FE face values are those from periodic boundary
                  fe_face_values.reinit(cell, face_id);
                  fe_face_values_periodic.reinit(periodic_cell, periodic_face);

                  for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                       ++f_q_point)
                    {
                      // Finding the normal vector of the boundary face
                      Tensor<1, dim> normal_vector =
                        fe_face_values.normal_vector(f_q_point);
                      Tensor<1, dim> periodic_normal_vector =
                        fe_face_values_periodic.normal_vector(f_q_point);

                      // Finding a point on the boundary face
                      Point<dim> quad_point =
                        fe_face_values.quadrature_point(0);
                      Point<dim> periodic_quad_point =
                        fe_face_values_periodic.quadrature_point(0);

                      // Storing these information for the current cell
                      boundary_cells_info_struct<dim> boundary_information;
                      boundary_information.cell = cell;
                      boundary_information.boundary_id =
                        cell->face(face_id)->boundary_id();
                      boundary_information.global_face_id =
                        cell->face_index(face_id);
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_face = quad_point;

                      // Storing these information for the periodic cell
                      boundary_cells_info_struct<dim>
                        periodic_boundary_information;
                      periodic_boundary_information.cell        = periodic_cell;
                      periodic_boundary_information.boundary_id = periodic_face;
                      periodic_boundary_information.global_face_id =
                        periodic_cell->face_index(periodic_face);
                      periodic_boundary_information.normal_vector =
                        periodic_normal_vector;
                      periodic_boundary_information.point_on_face =
                        periodic_quad_point;


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
                          if (cell->is_locally_owned())
                            {
                              periodic_boundary_cells_information.insert(
                                {cell->face_index(face_id),
                                 std::make_pair(
                                   boundary_information,
                                   periodic_boundary_information)});
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

template <int dim>
void
PeriodicBoundariesManipulator<dim>::execute_particle_displacement(
  const std::map<
    int,
    std::pair<boundary_cells_info_struct<dim>, boundary_cells_info_struct<dim>>>
    &                                    boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto boundary_cells_information_iterator =
         boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator)
    {
      auto boundary_cells_content =
        boundary_cells_information_iterator->second.first;
      auto cell = boundary_cells_content.cell;

      auto periodic_boundary_cells_content =
        boundary_cells_information_iterator->second.second;

      // Finding particles located in the corresponding cell
      // (boundary_cells_content.cell)
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      const bool particles_exist_in_main_cell = !particles_in_cell.empty();

      // If the main cell is not empty
      if (particles_exist_in_main_cell)
        {
          for (auto particle = particles_in_cell.begin();
               particle != particles_in_cell.end();
               ++particle)
            {
              Point<3> particle_position;
              double   distance_with_face;

              if constexpr (dim == 3)
                {
                  particle_position = particle->get_location();
                  distance_with_face =
                    scalar_product(particle_position -
                                     boundary_cells_content.point_on_face,
                                   boundary_cells_content.normal_vector);
                }

              if constexpr (dim == 2)
                {
                  particle_position = point_nd_to_3d(particle->get_location());
                  // distance_with_face =
                  // scalar_product(particle_position -
                  //               point_nd_to_3d(
                  //               boundary_cells_content.point_on_face),
                  //          boundary_cells_content.normal_vector);
                }



              if (distance_with_face >= 0.0)
                {
                  std::cout << "ID : " << particle->get_id() << std::endl;
                  std::cout << "distance with face : " << distance_with_face
                            << std::endl;
                  std::cout << "Old particle position " << particle_position
                            << std::endl;
                  std::cout << "reference position : "
                            << particle->get_reference_location() << std::endl;
                  Tensor<1, 3> distance_between_faces;

                  for (int d = 0; d < 3; ++d)
                    {
                      // Calculate distance between boundary faces
                      distance_between_faces[d] =
                        periodic_boundary_cells_content.point_on_face[d] -
                        boundary_cells_content.point_on_face[d];

                      particle_position[d] += distance_between_faces[d];
                    }

                  if constexpr (dim == 3)
                    particle->set_location(particle_position);

                  if constexpr (dim == 2)
                    {
                      Point<2> position_2d;
                      position_2d[0] = particle_position[0];
                      position_2d[1] = particle_position[1];
                      particle->set_location(position_2d);
                    }
                  std::cout << "New particle position " << particle_position
                            << std::endl;
                }
            }
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::find_cell_neighbors_and_periodic(
  const parallel::distributed::Triangulation<dim> &triangulation,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_local_neighbor_list,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_ghost_neighbor_list,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_local_periodic_neighbor_list,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_ghost_periodic_neighbor_list)
{
  // Find cell neighbors as no periodic boundaries
  find_cell_neighbors(triangulation,
                      cells_local_neighbor_list,
                      cells_ghost_neighbor_list);

  std::vector<typename Triangulation<dim>::active_cell_iterator>
    local_periodic_neighbor_vector;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    ghost_periodic_neighbor_vector;

  // This vector is used to avoid repetition of adjacent cells. For instance if
  // cell B is recognized as the neighbor of cell A, cell A will not be added to
  // the neighbor list of cell B again. This is done using the totall_cell_list
  // vector
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    totall_cell_list;

  // For each cell, the cell vertices are found and used to find adjacent cells.
  // The reason is to find the cells located on the corners of the main cell.
  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  unsigned int cell_number_iterator = 0;

  // Looping over cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
          // The first element of each vector is the cell itself.
          local_periodic_neighbor_vector.push_back(cell);
          totall_cell_list.push_back(cell);

          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              if (cell->face(face_id)->at_boundary() &&
                  cell->has_periodic_neighbor(face_id))
                {
                  for (unsigned int vertex = 0;
                       vertex < GeometryInfo<dim>::vertices_per_cell;
                       ++vertex)
                    {
                      for (const auto &periodic_neighbor :
                           v_to_c[cell->vertex_index(vertex)])
                        {
                          if (periodic_neighbor->is_locally_owned())
                            {
                              auto search_iterator =
                                std::find(totall_cell_list.begin(),
                                          totall_cell_list.end(),
                                          periodic_neighbor);
                              auto local_search_iterator = std::find(
                                local_periodic_neighbor_vector.begin(),
                                local_periodic_neighbor_vector.end(),
                                periodic_neighbor);

                              if (search_iterator == totall_cell_list.end() &&
                                  local_search_iterator ==
                                    local_periodic_neighbor_vector.end())
                                {
                                  local_periodic_neighbor_vector.push_back(
                                    periodic_neighbor);
                                }

                              // If the neighbor cell is a ghost, it should be
                              // added in the ghost_neighbor_vector container
                            }
                          else if (periodic_neighbor->is_ghost())
                            {
                              auto ghost_search_iterator = std::find(
                                ghost_periodic_neighbor_vector.begin(),
                                ghost_periodic_neighbor_vector.end(),
                                periodic_neighbor);
                              if (ghost_search_iterator ==
                                  ghost_periodic_neighbor_vector.end())
                                {
                                  if (ghost_periodic_neighbor_vector.empty())
                                    {
                                      ghost_periodic_neighbor_vector.push_back(
                                        cell);
                                    }

                                  ghost_periodic_neighbor_vector.push_back(
                                    periodic_neighbor);
                                }
                            }
                        }
                    }
                }
            }
        }
      if (!local_periodic_neighbor_vector.empty())
        cells_local_periodic_neighbor_list.push_back(
          local_periodic_neighbor_vector);
      if (!ghost_periodic_neighbor_vector.empty())
        cells_ghost_periodic_neighbor_list.push_back(
          ghost_periodic_neighbor_vector);
      local_periodic_neighbor_vector.clear();
      ghost_periodic_neighbor_vector.clear();
      ++cell_number_iterator;
    }
}


template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;