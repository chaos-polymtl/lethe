#include <dem/particle_wall_fine_search.h>

#include <boost/math/special_functions/sign.hpp>


using namespace dealii;

template <int dim>
ParticleWallFineSearch<dim>::ParticleWallFineSearch()
{}

template <int dim>
void ParticleWallFineSearch<dim>::particle_wall_fine_search(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  types::boundary_id,
                                  unsigned int>>>
    &particle_wall_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact)
{
  // Iterating over contact candidates from broad search and adding the pairs to
  // the particle_wall_pairs_in_contact
  for (auto const &[particle_id, particle_pair_candidates] :
       particle_wall_contact_pair_candidates)
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

          // Setting tangential overlap of the new particle-wall contact
          // pair equal to zero
          Tensor<1, 3> tangential_overlap;
          tangential_overlap[0] = 0.0;
          tangential_overlap[1] = 0.0;
          tangential_overlap[2] = 0.0;


          Tensor<1, 3> normal_vector_3d;
          if constexpr (dim == 3)
            normal_vector_3d = normal_vector;

          if constexpr (dim == 2)
            {
              normal_vector_3d[0] = normal_vector[0];
              normal_vector_3d[1] = normal_vector[1];
              normal_vector_3d[2] = 0.0;
            }

          Point<3> point_on_boundary_3d;
          if constexpr (dim == 3)
            point_on_boundary_3d = point_on_boundary;

          if constexpr (dim == 2)
            {
              point_on_boundary_3d[0] = point_on_boundary[0];
              point_on_boundary_3d[1] = point_on_boundary[1];
              point_on_boundary_3d[2] = 0.0;
            }


          // Adding contact info to the sample to
          // particle_wall_contact_info_struct
          particle_wall_contact_info_struct<dim> contact_info;
          contact_info.particle                 = particle;
          contact_info.normal_vector            = normal_vector_3d;
          contact_info.normal_overlap           = .0;
          contact_info.normal_relative_velocity = .0;
          contact_info.point_on_boundary        = point_on_boundary_3d;
          contact_info.boundary_id =
            std::get<3>(particle_pair_candidate_content);
          contact_info.tangential_overlap           = tangential_overlap;
          contact_info.tangential_relative_velocity = .0;
          contact_info.global_face_id =
            std::get<4>(particle_pair_candidate_content);

          particle_wall_pairs_in_contact[particle_id].insert(
            {face_id, contact_info});
        }
    }
}

template <int dim>
void
ParticleWallFineSearch<dim>::particle_floating_wall_fine_search(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>>
    &                                               pfw_contact_candidates,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double &                                    simulation_time,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &pfw_pairs_in_contact)
{
  // Reading floating wall properties
  std::vector<Point<dim>> point_on_wall =
    floating_wall_properties.points_on_walls;
  std::vector<Tensor<1, dim>> wall_normal_vector =
    floating_wall_properties.floating_walls_normal_vectors;

  // Iterating over contact candidates from broad search and adding the pairs to
  // the pfw_pairs_in_contact
  for (auto const &[particle_id, particle_pair_candidates] :
       pfw_contact_candidates)
    {
      for (auto particle_pair_candidate_iterator =
             particle_pair_candidates.begin();
           particle_pair_candidate_iterator != particle_pair_candidates.end();
           ++particle_pair_candidate_iterator)
        {
          // Getting the floating wall id once to improve efficiency
          unsigned int floating_wall_id =
            particle_pair_candidate_iterator->first;

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
                {
                  point_on_floating_wall_3d[0] = point_on_floating_wall[0];
                  point_on_floating_wall_3d[1] = point_on_floating_wall[1];
                  point_on_floating_wall_3d[2] = 0.0;
                }

              // Check to see on which side of the wall the particle is located:

              // Finding connecting vector from defined point on the boundary
              // wall to the particle location
              Tensor<1, dim> connecting_vector =
                particle->get_location() - point_on_floating_wall;
              int inner_product_sign =
                boost::math::sign(connecting_vector * normal_vector);

              // If the cell is located on the opposite side of the defined
              // normal vector, the normal vector of the cell should be reversed
              if (inner_product_sign < 0)
                {
                  normal_vector = -1 * normal_vector;
                }

              Tensor<1, 3> normal_vector_3d;
              if constexpr (dim == 3)
                normal_vector_3d = normal_vector;

              if constexpr (dim == 2)
                {
                  normal_vector_3d[0] = normal_vector[0];
                  normal_vector_3d[1] = normal_vector[1];
                  normal_vector_3d[2] = 0.0;
                }

              // Setting tangential overlap of the new particle-floating wall
              // contact pair equal to zero
              // VECTORIZE ***********
              Tensor<1, 3> tangential_overlap;
              tangential_overlap[0] = 0.0;
              tangential_overlap[1] = 0.0;
              tangential_overlap[2] = 0.0;


              // Creating a sample from the particle_wall_contact_info_struct
              // and adding contact info to the sample
              particle_wall_contact_info_struct<dim> contact_info;
              contact_info.particle                 = particle;
              contact_info.normal_vector            = normal_vector_3d;
              contact_info.normal_overlap           = .0;
              contact_info.normal_relative_velocity = .0;
              contact_info.point_on_boundary        = point_on_floating_wall_3d;
              // The boundary ID of floating walls is set to 100, it should be
              // modified after adding motion of floating walls
              contact_info.boundary_id                  = 100;
              contact_info.global_face_id               = 0;
              contact_info.tangential_overlap           = tangential_overlap;
              contact_info.tangential_relative_velocity = .0;


              pfw_pairs_in_contact[particle_id].insert(
                {floating_wall_id, contact_info});
            }
        }
    }
}

template class ParticleWallFineSearch<2>;
template class ParticleWallFineSearch<3>;
