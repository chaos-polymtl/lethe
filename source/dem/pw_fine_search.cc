#include <dem/pw_fine_search.h>

using namespace dealii;

template <int dim>
PWFineSearch<dim>::PWFineSearch()
{}

template <int dim>
void PWFineSearch<dim>::particle_wall_fine_search(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  types::boundary_id,
                                  unsigned int>>> &pw_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    &pw_pairs_in_contact)
{
  // Iterating over contact candidates from broad search and adding the pairs to
  // the pw_pairs_in_contact
  for (auto const &[particle_id, particle_pair_candidates] :
       pw_contact_pair_candidates)
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
          Tensor<1, dim> tangential_overlap;
          tangential_overlap[0] = 0.0;
          tangential_overlap[1] = 0.0;
          if (dim == 3)
            {
              tangential_overlap[2] = 0.0;
            }

          // Adding contact info to the sample to pw_contact_info_struct
          pw_contact_info_struct<dim> contact_info;
          contact_info.particle                 = particle;
          contact_info.normal_vector            = normal_vector;
          contact_info.normal_overlap           = .0;
          contact_info.normal_relative_velocity = .0;
          contact_info.point_on_boundary        = point_on_boundary;
          contact_info.boundary_id =
            std::get<3>(particle_pair_candidate_content);
          contact_info.tangential_overlap           = tangential_overlap;
          contact_info.tangential_relative_velocity = .0;
          contact_info.global_face_id =
            std::get<4>(particle_pair_candidate_content);

          pw_pairs_in_contact[particle_id].insert({face_id, contact_info});
        }
    }
}

template <int dim>
void
PWFineSearch<dim>::particle_floating_wall_fine_search(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>>
    &                                               pfw_contact_candidates,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double &                                    simulation_time,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
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

              // Setting tangential overlap of the new particle-floating wall
              // contact pair equal to zero
              Tensor<1, dim> tangential_overlap;
              tangential_overlap[0] = 0.0;
              tangential_overlap[1] = 0.0;
              if (dim == 3)
                {
                  tangential_overlap[2] = 0.0;
                }

              // Creating a sample from the pw_contact_info_struct and adding
              // contact info to the sample
              pw_contact_info_struct<dim> contact_info;
              contact_info.particle                 = particle;
              contact_info.normal_vector            = normal_vector;
              contact_info.normal_overlap           = .0;
              contact_info.normal_relative_velocity = .0;
              contact_info.point_on_boundary        = point_on_floating_wall;
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

template class PWFineSearch<2>;
template class PWFineSearch<3>;
