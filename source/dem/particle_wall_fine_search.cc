#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_wall_fine_search.h>

#include <boost/math/special_functions/sign.hpp>


using namespace dealii;

template <int dim>
ParticleWallFineSearch<dim>::ParticleWallFineSearch()
{}

template <int dim>
void
ParticleWallFineSearch<dim>::particle_wall_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_pair_candidates,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_wall_pairs_in_contact)
{
  // Iterating over contact candidates from broad search and adding the pairs to
  // the particle_wall_pairs_in_contact
  for (auto const &[particle_id, particle_pair_candidates] :
       particle_wall_contact_pair_candidates)
    {
      if (!particle_pair_candidates.empty())
        {
          for (auto const &[face_id, particle_pair_candidate_content] :
               particle_pair_candidates)
            {
              auto particle = std::get<0>(particle_pair_candidate_content);

              // Normal vector of the boundary and a point on the boudary are
              // defined as local parameters
              auto normal_vector = std::get<1>(particle_pair_candidate_content);
              Point<dim> point_on_boundary =
                std::get<2>(particle_pair_candidate_content);

              Tensor<1, 3> normal_vector_3d;
              if constexpr (dim == 3)
                normal_vector_3d = normal_vector;

              if constexpr (dim == 2)
                normal_vector_3d = tensor_nd_to_3d(normal_vector);

              Point<3> point_on_boundary_3d;
              if constexpr (dim == 3)
                point_on_boundary_3d = point_on_boundary;

              if constexpr (dim == 2)
                point_on_boundary_3d = point_nd_to_3d(point_on_boundary);

              particle_wall_pairs_in_contact[particle_id].emplace(
                face_id,
                particle_wall_contact_info<dim>(
                  particle,
                  normal_vector_3d,
                  point_on_boundary_3d,
                  std::get<3>(particle_pair_candidate_content),
                  std::get<4>(particle_pair_candidate_content)));
            }
        }
    }
}

template <int dim>
void
ParticleWallFineSearch<dim>::particle_floating_wall_fine_search(
  const typename DEM::dem_data_structures<
    dim>::particle_floating_wall_candidates &particle_floating_wall_candidates,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_floating_wall_pairs_in_contact)
{
  // Reading floating wall properties
  std::vector<Point<dim>> point_on_wall =
    floating_wall_properties.points_on_walls;
  std::vector<Tensor<1, dim>> wall_normal_vector =
    floating_wall_properties.floating_walls_normal_vectors;

  // Iterating over contact candidates from broad search and adding the pairs to
  // the particle_floating_wall_pairs_in_contact
  for (auto const &[particle_id, particle_pair_candidates] :
       particle_floating_wall_candidates)
    {
      if (!particle_pair_candidates.empty())
        {
          for (auto particle_pair_candidate_iterator =
                 particle_pair_candidates.begin();
               particle_pair_candidate_iterator !=
               particle_pair_candidates.end();
               ++particle_pair_candidate_iterator)
            {
              // Getting the floating wall id once to improve efficiency
              auto floating_wall_id = particle_pair_candidate_iterator->first;

              // Checking simulation time for temporary floating walls
              if (simulation_time >=
                    floating_wall_properties.time_start[floating_wall_id] &&
                  simulation_time <=
                    floating_wall_properties.time_end[floating_wall_id])
                {
                  // Reading particle, normal vector and point on wall once to
                  // improve efficiency
                  auto particle = particle_pair_candidate_iterator->second;
                  Tensor<1, dim> normal_vector =
                    wall_normal_vector[floating_wall_id];
                  Point<dim> point_on_floating_wall =
                    point_on_wall[floating_wall_id];

                  Point<3> point_on_floating_wall_3d;
                  if constexpr (dim == 3)
                    point_on_floating_wall_3d = point_on_floating_wall;

                  if constexpr (dim == 2)
                    point_on_floating_wall_3d =
                      point_nd_to_3d(point_on_floating_wall);

                  // Check to see on which side of the wall the particle is
                  // located:

                  // Finding connecting vector from defined point on the
                  // boundary wall to the particle location
                  Tensor<1, dim> connecting_vector =
                    particle->get_location() - point_on_floating_wall;
                  int inner_product_sign =
                    boost::math::sign(connecting_vector * normal_vector);

                  // If the cell is located on the opposite side of the defined
                  // normal vector, the normal vector of the cell should be
                  // reversed
                  if (inner_product_sign < 0)
                    {
                      normal_vector = -1 * normal_vector;
                    }

                  Tensor<1, 3> normal_vector_3d;
                  if constexpr (dim == 3)
                    normal_vector_3d = normal_vector;

                  if constexpr (dim == 2)
                    normal_vector_3d = tensor_nd_to_3d(normal_vector);

                  particle_floating_wall_pairs_in_contact[particle_id].emplace(
                    floating_wall_id,
                    particle_wall_contact_info<dim>(particle,
                                                    normal_vector_3d,
                                                    point_on_floating_wall_3d,
                                                    100,
                                                    0));
                }
            }
        }
    }
}


template <int dim>
void
ParticleWallFineSearch<dim>::particle_floating_mesh_fine_search(
  const typename DEM::dem_data_structures<
    dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
    &particle_floating_mesh_in_contact)
{
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_contact_candidates.size();
       ++solid_counter)
    {
      auto &particle_floating_mesh_element =
        particle_floating_mesh_in_contact[solid_counter];

      auto &candidates =
        particle_floating_mesh_contact_candidates[solid_counter];

      for (auto const &[cut_cell_key, candidate_particles] : candidates)
        {
          if (!candidate_particles.empty())
            {
              for (auto &particle_floating_mesh_candidate_iterator :
                   candidate_particles)
                {
                  particle_floating_mesh_element[cut_cell_key].emplace(
                    particle_floating_mesh_candidate_iterator.first,
                    particle_wall_contact_info<dim>(
                      particle_floating_mesh_candidate_iterator.second));
                }
            }
        }
    }
}

template class ParticleWallFineSearch<2>;
template class ParticleWallFineSearch<3>;
