#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/periodic_boundaries_manupulator.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

using namespace dealii;

template <int dim>
PeriodicBoundariesManipulator<dim>::PeriodicBoundariesManipulator()
{}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::get_boundary_info(
  typename Triangulation<dim>::cell_iterator cell,
  unsigned int                               face_id,
  boundary_cells_info_struct<dim> &          boundary_information)
{
  // Initialize a simple quadrature for on the system. This will be used to
  // obtain a single sample point on the boundary faces
  const FE_Q<dim>   fe(1);
  QGauss<dim - 1>   face_quadrature_formula(1);
  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);

  // Store information
  boundary_information.cell           = cell;
  boundary_information.boundary_id    = cell->face(face_id)->boundary_id();
  boundary_information.global_face_id = cell->face_index(face_id);

  // Finding the normal vector of the boundary face and point
  fe_face_values.reinit(cell, face_id);
  Tensor<1, dim> normal_vector = fe_face_values.normal_vector(0);
  Point<dim>     quad_point    = fe_face_values.quadrature_point(0);

  boundary_information.normal_vector = normal_vector;
  boundary_information.point_on_face = quad_point;
}



template <int dim>
void
PeriodicBoundariesManipulator<dim>::map_periodic_cells(
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


  /*
  // Get map of periodic cells at periodic boundaries as :
  // map[<cell_iterator, cell_face_id>,
  //     <<periodic_cell_iterator, periodic_cell_face_id>, orientation>]
  std::map<std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>,
           std::pair<std::pair<typename Triangulation<dim>::cell_iterator,
                               unsigned int>,
                     std::bitset<3>>>
    periodic_face_map = triangulation.get_periodic_face_map(); */



  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // Iterating over the faces of each cell
      for (int face_id = 0; face_id < int(GeometryInfo<dim>::faces_per_cell);
           ++face_id)
        {
          if (cell->has_periodic_neighbor(face_id))
            {
              // Check to have unique pair
              if (!periodic_cell_pair.count(
                    cell->periodic_neighbor(face_id)->active_cell_index()))
                {
                  // Save boundaries information related to cell, face and
                  // boundary ids Cell on boundary tagged as outlet
                  boundary_cells_info_struct<dim> boundary_information;

                  get_boundary_info(cell, face_id, boundary_information);

                  boundary_cells_info_struct<dim> periodic_boundary_information;
                  typename Triangulation<dim>::active_cell_iterator
                    periodic_cell = cell->periodic_neighbor(face_id);

                  get_boundary_info(periodic_cell,
                                    cell->periodic_neighbor_face_no(face_id),
                                    periodic_boundary_information);

                  auto information_search_element =
                    std::make_pair(face_id, cell);
                  auto search_iterator = std::find(search_vector.begin(),
                                                   search_vector.end(),
                                                   information_search_element);
                  if (search_iterator == search_vector.end())
                    {
                      if (cell->is_locally_owned())
                        {
                          search_vector.push_back(information_search_element);
                          periodic_boundary_cells_information.insert(
                            {boundary_information.cell->active_cell_index(),
                             std::make_pair(boundary_information,
                                            periodic_boundary_information)});
                          periodic_cell_pair.insert(
                            {boundary_information.cell->active_cell_index(),
                             periodic_boundary_information.cell
                               ->active_cell_index()});

                          /*std::cout << "--" <<
                          boundary_information.cell->active_cell_index()
                                    << std::endl;
                          std::cout << "-" <<
                          periodic_cell_pair[boundary_information.cell->active_cell_index()]
                                    << std::endl;
                          std::cout << boundary_information.global_face_id <<
                          std::endl
                            ;
                          std::cout << boundary_information.point_on_face <<
                          std::endl; std::cout <<
                          periodic_boundary_information.global_face_id
                                    << std::endl;
                          std::cout <<
                          periodic_boundary_information.point_on_face
                                    << std::endl;*/
                        }
                    }
                }
            }
        }


      unsigned int first_cell_id = 0;
      if (cell->active_cell_index() == first_cell_id &&
          cell->is_locally_owned())
        {
          translation_value = periodic_boundary_cells_information[first_cell_id]
                                .first.point_on_face[direction] -
                              periodic_boundary_cells_information[first_cell_id]
                                .second.point_on_face[direction];
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::move_particles(
  boundary_cells_info_struct<dim> &cell_1,
  boundary_cells_info_struct<dim> &cell_2,
  typename Particles::ParticleHandler<dim>::particle_iterator_range
    &particles_in_cell)
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
          // Calculate distance between particle position and the cell
          distance_with_face =
            scalar_product(particle_position - cell_1.point_on_face,
                           cell_1.normal_vector);
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


      // Check if distance is positive, if so, it is toward the norm,
      // which means outside of cell.
      // If particle is outside of the domain of cell on the periodic
      // boundary, it is move in the periodic cell.
      if (distance_with_face >= 0.0)
        {
          std::cout << "ID : " << particle->get_id() << std::endl;
          std::cout << "distance with face : " << distance_with_face
                    << std::endl;
          std::cout << "Old particle position " << particle_position
                    << std::endl;
          std::cout << "reference position : "
                    << particle->get_reference_location() << std::endl;

          double distance_between_faces = cell_2.point_on_face[direction] -
                                          cell_1.point_on_face[direction];

          std::cout << distance_between_faces << std::endl;

          // Move particle outside the current cell to inside the
          // periodic cell.
          particle_position[direction] += distance_between_faces;

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

template <int dim>
void
PeriodicBoundariesManipulator<dim>::execute_particle_displacement(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto boundary_cells_information_iterator =
         periodic_boundary_cells_information.begin();
       boundary_cells_information_iterator !=
       periodic_boundary_cells_information.end();
       ++boundary_cells_information_iterator)
    {
      auto boundary_cells_content =
        boundary_cells_information_iterator->second.first;
      auto cell = boundary_cells_content.cell;
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      auto periodic_boundary_cells_content =
        boundary_cells_information_iterator->second.second;
      auto periodic_cell = periodic_boundary_cells_content.cell;
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_periodic_cell =
          particle_handler.particles_in_cell(periodic_cell);


      const bool particles_exist_in_main_cell = !particles_in_cell.empty();
      const bool particles_exist_in_periodic_cell =
        !particles_in_periodic_cell.empty();

      // If the main cell is not empty
      if (particles_exist_in_main_cell)
        {
          move_particles(boundary_cells_content,
                         periodic_boundary_cells_content,
                         particles_in_cell);
        }
      if (particles_exist_in_periodic_cell)
        {
          move_particles(periodic_boundary_cells_content,
                         boundary_cells_content,
                         particles_in_periodic_cell);
        }
    }
}



/*
template <int dim>
void
PeriodicBoundariesManipulator<dim>::find_periodic_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_local_periodic_neighbor_list,
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    &cells_ghost_periodic_neighbor_list)
{
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    local_periodic_neighbor_vector;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    ghost_periodic_neighbor_vector;

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
          if(cell->active_cell_index())

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

template <int dim>
void
PeriodicBoundariesManipulator<dim>::find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const std::map<
    int,
    std::pair<boundary_cells_info_struct<dim>, boundary_cells_info_struct<dim>>>
    &periodic_boundaries_info,
  const std::vector<
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
    *cells_local_periodic_neighbor_list,
  const std::vector<
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
    *cells_ghost_periodic_neighbor_list,
  std::unordered_map<types::particle_index,
                     std::vector<std::pair<types::particle_index, Point<3>>>>
    &local_periodic_contact_pair_candidates,
  std::unordered_map<types::particle_index,
                     std::vector<std::pair<types::particle_index, Point<3>>>>
    &ghost_periodic_contact_pair_candidates)
{
  local_periodic_contact_pair_candidates.clear();

  for (auto cell_periodic_neighbor_list_iterator =
         cells_local_periodic_neighbor_list->begin();
       cell_periodic_neighbor_list_iterator !=
       cells_local_periodic_neighbor_list->end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      // Check to see if the main cell has any particles
      if (particles_exist_in_main_cell)
        {
          // Going through periodic neighbor cells of the main cell
          ++cell_periodic_neighbor_iterator;

          for (; cell_periodic_neighbor_iterator !=
                 cell_periodic_neighbor_list_iterator->end();
               ++cell_periodic_neighbor_iterator)
            {
              // Defining iterator on local particles in the neighbor cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_neighbor_cell =
                  particle_handler.particles_in_cell(
                    *cell_periodic_neighbor_iterator);

              // Capturing particle pairs, the first particle in the main
              // cell and the second particle in the neighbor cells
              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_main_cell_iterator =
                       particles_in_main_cell.begin();
                   particles_in_main_cell_iterator !=
                   particles_in_main_cell.end();
                   ++particles_in_main_cell_iterator)
                {
                  std::vector<std::pair<types::particle_index, Point<3>>>
                    &particle_periodic_candidate_container =
                      local_periodic_contact_pair_candidates
                        [particles_in_main_cell_iterator->get_id()];
                  if (particle_periodic_candidate_container.empty())
                    particle_periodic_candidate_container.reserve(40);

                  for (typename Particles::ParticleHandler<
                         dim>::particle_iterator_range::iterator
                         particles_in_periodic_neighbor_cell_iterator =
                           particles_in_periodic_neighbor_cell.begin();
                       particles_in_periodic_neighbor_cell_iterator !=
                       particles_in_periodic_neighbor_cell.end();
                       ++particles_in_periodic_neighbor_cell_iterator)
                    {
                      // Get the virtual periodic location of particle for main
                      // particle
                      Point<3> periodic_location(
                        particles_in_periodic_neighbor_cell_iterator
                          ->get_location());

                      for (int face_id = 0;
                           face_id < int(GeometryInfo<dim>::faces_per_cell);
                           ++face_id)
                        {
                          if (cell_periodic_neighbor_iterator->face(face_id)
                                ->at_boundary() &&
                              cell_periodic_neighbor_iterator
                                ->has_periodic_neighbor(face_id))
                            {
                            }



                          particle_periodic_candidate_container.emplace_back(
                            particles_in_periodic_neighbor_cell_iterator
                              ->get_id());
                        }
                    }
                }
            }
        }
    }
}
}*/


template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;