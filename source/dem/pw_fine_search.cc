#include <dem/pw_fine_search.h>

using namespace dealii;

template <int dim>
PWFineSearch<dim>::PWFineSearch()
{}

template <int dim>
void
PWFineSearch<dim>::pw_Fine_Search(
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    &pw_contact_pair_candidates,
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    &pw_pairs_in_contact)
{
  // Now iterating over contact candidates from broad search. If a particle-wall
  // pair is in contact (distance > 0) and does not exist in the
  // pw_pairs_in_contact, it is added to the pw_pairs_in_contact
  for (auto map_iterator = pw_contact_pair_candidates.begin();
       map_iterator != pw_contact_pair_candidates.end();
       ++map_iterator)
    {
      auto particle_id              = map_iterator->first;
      auto particle_pair_candidates = &map_iterator->second;

      for (auto particle_pair_candidate_iterator =
             particle_pair_candidates->begin();
           particle_pair_candidate_iterator != particle_pair_candidates->end();
           ++particle_pair_candidate_iterator)
        {
          // Get the particle and face id from the vector and the total array
          // view to the particle properties once to improve efficiency
          int  face_id = particle_pair_candidate_iterator->first;
          auto particle_pair_candidate_content =
            particle_pair_candidate_iterator->second;
          auto particle = std::get<0>(particle_pair_candidate_content);

          // Normal vector of the boundary and a point on the boudary are
          // defined as local parameters
          Tensor<1, dim> normal_vector =
            std::get<1>(particle_pair_candidate_content);
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

          // Creating a sample from the pw_contact_info_struct and adding
          // contact info to the sample
          pw_contact_info_struct<dim> contact_info;
          contact_info.particle           = particle;
          contact_info.normal_vector      = normal_vector;
          contact_info.point_on_boundary  = point_on_boundary;
          contact_info.tangential_overlap = tangential_overlap;

          pw_pairs_in_contact[particle_id].insert({face_id, contact_info});
        }
    }
}

template class PWFineSearch<2>;
template class PWFineSearch<3>;
