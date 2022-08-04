#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/periodic_boundaries_manipulator.h>

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
  periodic_boundary_cells_information.clear();
  global_periodic_cell_pair.clear();

  // Iterating over the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->at_boundary())
        {
          // Iterating over the faces of each cell
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              if (cell->has_periodic_neighbor(face_id))
                {
                  // Check the global cell index key prior having unique pair
                  if (!global_periodic_cell_pair.count(
                        cell->periodic_neighbor(face_id)->active_cell_index()))
                    {
                      // Save boundaries information related to the cell on
                      // boundary tagged as outlet
                      boundary_cells_info_struct<dim> boundary_information;
                      get_boundary_info(cell, face_id, boundary_information);

                      // Save boundaries information related to the cell on
                      // boundary tagged as periodic
                      boundary_cells_info_struct<dim>
                        periodic_boundary_information;
                      typename Triangulation<dim>::active_cell_iterator
                        periodic_cell = cell->periodic_neighbor(face_id);
                      get_boundary_info(periodic_cell,
                                        cell->periodic_neighbor_face_no(
                                          face_id),
                                        periodic_boundary_information);

                      // Store both cell information in map with cell id at
                      // outlet as key
                      if (cell->is_locally_owned())
                        {
                          periodic_boundary_cells_information.insert(
                            {boundary_information.cell->active_cell_index(),
                             std::make_pair(boundary_information,
                                            periodic_boundary_information)});
                        }

                      // Map of the periodic cell related to the cell at outlet
                      global_periodic_cell_pair.insert(
                        {boundary_information.cell->active_cell_index(),
                         periodic_boundary_information.cell
                           ->active_cell_index()});

                    }
                }
            }
        }
    }

  std::cout << " is empty " << periodic_boundary_cells_information.empty()
            << std::endl;
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::check_and_move_particles(
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
      double   distance_with_face = 0;

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
          double distance_between_faces =
            cell_2.point_on_face[direction] - cell_1.point_on_face[direction];

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
        }
    }
}

template <int dim>
void
PeriodicBoundariesManipulator<dim>::execute_particle_displacement(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  if (!periodic_boundary_cells_information.empty())
    {
      for (auto boundary_cells_information_iterator =
             periodic_boundary_cells_information.begin();
           boundary_cells_information_iterator !=
           periodic_boundary_cells_information.end();
           ++boundary_cells_information_iterator)
        {
          // Get the cell and periodic cell content from map
          auto boundary_cells_content =
            boundary_cells_information_iterator->second.first;
          auto cell = boundary_cells_content.cell;

          auto periodic_boundary_cells_content =
            boundary_cells_information_iterator->second.second;
          auto periodic_cell = periodic_boundary_cells_content.cell;

          /*std::cout << "Cell index key: "
                    << boundary_cells_information_iterator->first << std::endl;
          std::cout << "Cell ID:        " << cell->active_cell_index()
                    << std::endl;
          std::cout << "Periodic cell ID: "
                    << periodic_cell->active_cell_index() << std::endl;
          std::cout << "i:" << i << std::endl;*/

          if (cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                         particles_in_cell = particle_handler.particles_in_cell(cell);
              const bool particles_exist_in_main_cell =
                !particles_in_cell.empty();

              if (particles_exist_in_main_cell)
                {
                  check_and_move_particles(boundary_cells_content,
                                           periodic_boundary_cells_content,
                                           particles_in_cell);
                }
            }

          if (periodic_cell->is_locally_owned())
            {
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_cell =
                  particle_handler.particles_in_cell(periodic_cell);

              const bool particles_exist_in_periodic_cell =
                !particles_in_periodic_cell.empty();

              if (particles_exist_in_periodic_cell)
                {
                  check_and_move_particles(periodic_boundary_cells_content,
                                           boundary_cells_content,
                                           particles_in_periodic_cell);
                }
            }
        }
    }
}

template class PeriodicBoundariesManipulator<2>;
template class PeriodicBoundariesManipulator<3>;