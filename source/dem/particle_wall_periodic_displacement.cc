#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_wall_periodic_displacement.h>

using namespace dealii;

template <int dim>
ParticleWallPeriodicDisplacement<dim>::ParticleWallPeriodicDisplacement()
{}

template <int dim>
void
ParticleWallPeriodicDisplacement<dim>::execute_particle_displacement(
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
                  //distance_with_face =
                    //scalar_product(particle_position -
                      //               point_nd_to_3d(
                        //               boundary_cells_content.point_on_face),
                         //          boundary_cells_content.normal_vector);
                }


              std::cout << "ID : " << particle->get_id() << std::endl;
              std::cout << "distance with face : " << distance_with_face << std::endl;

              if (distance_with_face >= 0.0)
                {
                  std::cout << "Old particle position " << particle_position << std::endl;
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
                  std::cout << "New particle position " << particle_position << std::endl;
                }
            }
        }
    }
}

template class ParticleWallPeriodicDisplacement<2>;
template class ParticleWallPeriodicDisplacement<3>;