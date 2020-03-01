#include <dem/pp_broad_search.h>

using namespace dealii;

template <int dim> PPBroadSearch<dim>::PPBroadSearch() {}

template <int dim>
std::vector<std::pair<Particles::ParticleIterator<dim>,
                      Particles::ParticleIterator<dim>>>
PPBroadSearch<dim>::find_PP_Contact_Pairs(
    Particles::ParticleHandler<dim> &particle_handler,
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
        &cell_Neighbor_List) {
  // A vector of pairs which contains all the canditates for particle-particle
  // collision at each time step
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
      contact_pair_candidates;

  // A counter for counting the cells
  int cell_number_counter = 0;

  // Looping over cell_neighbor_list
  for (auto cell_neighbor_list_iterator = cell_Neighbor_List.begin();
       cell_neighbor_list_iterator != cell_Neighbor_List.end();
       ++cell_neighbor_list_iterator, ++cell_number_counter) {

    // The main cell
    auto neighbor_cell_iterator =
        cell_Neighbor_List[cell_number_counter].begin();

    // Particles in the main cell
    typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
            particle_handler.particles_in_cell(*neighbor_cell_iterator);

    // finding collision pairs in the main cell, particle counter starts from 1,
    // becasue each particle will not be considered as collision partner with
    // itself
    int particle_in_main_cell_counter = 1;
    for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_main_cell_iterator =
                 particles_in_main_cell.begin();
         particles_in_main_cell_iterator != particles_in_main_cell.end();
         ++particles_in_main_cell_iterator, ++particle_in_main_cell_counter) {

      typename Particles::ParticleHandler<dim>::particle_iterator_range::
          iterator particles_in_main_cell_iterator_two =
              particles_in_main_cell.begin();

      // Advancing the second iterator to capture all the particle pairs in the
      // main cell
      std::advance(particles_in_main_cell_iterator_two,
                   particle_in_main_cell_counter);

      // while (particles_in_main_cell_iterator_two !=
      //  particles_in_main_cell.end()) {
      for (particles_in_main_cell_iterator_two;
           particles_in_main_cell_iterator_two != particles_in_main_cell.end();
           ++particles_in_main_cell_iterator_two) {
        // Capturing all the particle pairs in the main cell
        auto contact_pair = std::make_pair(particles_in_main_cell_iterator,
                                           particles_in_main_cell_iterator_two);
        contact_pair_candidates.push_back(contact_pair);
        //++particles_in_main_cell_iterator_two;
      }
    }

    // Going through neighbor cells of the main cell
    ++neighbor_cell_iterator;

    // while (neighbor_cell_iterator !=
    //   cell_Neighbor_List[cell_number_counter].end()) {
    //************* unused warning check
    for (neighbor_cell_iterator; neighbor_cell_iterator !=
                                 cell_Neighbor_List[cell_number_counter].end();
         ++neighbor_cell_iterator) {
      // Defining iterator on particles in the neighbor cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
          particles_in_neighbor_cell =
              particle_handler.particles_in_cell(*neighbor_cell_iterator);

      // Capturing particle pairs, the first particle in the main cell and the
      // second particle in the neighbor cells
      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
               iterator particles_in_main_cell_iterator =
                   particles_in_main_cell.begin();
           particles_in_main_cell_iterator != particles_in_main_cell.end();
           ++particles_in_main_cell_iterator) {
        for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
                 iterator particles_in_neighbor_cell_iterator =
                     particles_in_neighbor_cell.begin();
             particles_in_neighbor_cell_iterator !=
             particles_in_neighbor_cell.end();
             ++particles_in_neighbor_cell_iterator) {
          auto contact_pair =
              std::make_pair(particles_in_main_cell_iterator,
                             particles_in_neighbor_cell_iterator);
          contact_pair_candidates.push_back(contact_pair);
        }
      }

      //++neighbor_cell_iterator;
    }
  }
  return contact_pair_candidates;
}

// These are old EXTREMELY unefficient methods for broad search I had written
// earlier
// 2nd method:
/*
template <int dim>
std::vector<std::pair<Particles::ParticleIterator<dim>,
                      Particles::ParticleIterator<dim>>>
ContactSearch<dim>::findContactPairs(
  Particles::ParticleHandler<dim> &particle_handler,
  const Triangulation<dim> &       tr,
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    cellNeighborList)


{
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
      contactPairs;
  int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell, ++index)
    {
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particle_range = particle_handler.particles_in_cell(cell);


      for (auto cellIt = cellNeighborList[index].begin();
           cellIt != cellNeighborList[index].end();
           cellIt++)
        {
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particle_range2 = particle_handler.particles_in_cell(*cellIt);

          for (typename Particles::ParticleHandler<dim>::
                 particle_iterator_range::iterator partIter =
                   particle_range.begin();
               partIter != particle_range.end();
               ++partIter)
            {
              for (typename Particles::ParticleHandler<dim>::
                     particle_iterator_range::iterator partIter2 =
                       particle_range2.begin();
                   partIter2 != particle_range2.end();
                   ++partIter2)
                {
                  auto cPair  = std::make_pair(partIter2, partIter);
                  auto cPair2 = std::make_pair(partIter, partIter2);
                  auto it2 =
                    std::find(contactPairs.begin(), contactPairs.end(), cPair);
                  auto it3 =    | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |         2 |      19.7s |        22% |
| Solve                           |         2 |      3.03s |       3.4% |
| Setup dof system                |         1 |      3.97s |       4.5% |
+---------------------------------+-----------+------------+------------+
                    std::find(contactPairs.begin(), contactPairs.end(), cPair2);
                  if (it2 == contactPairs.end())
                    if (it3 == contactPairs.end())
                      if (partIter2 != partIter)
                        contactPairs.push_back(cPair);
                }
            }
        }
    }

  return contactPairs;
}
*/

// 1st method:
/*
{    | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble                        |         2 |      19.7s |        22% |
| Solve                           |         2 |      3.03s |       3.4% |
| Setup dof system                |         1 |      3.97s |       4.5% |
+---------------------------------+-----------+------------+------------+
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
    contactPairs;
  if (!contactPairs.empty())
    {
      contactPairs.clear();
    }
  typename Triangulation<dim>::active_cell_iterator currrentCell;
  for (auto particleIter = particle_handler.begin();
       particleIter != particle_handler.end();
       ++particleIter)
    {
      currrentCell =
        GridTools::find_active_cell_around_point(tr,
                                                 particleIter->get_location());
      auto it1 =
        std::find(totallCellList.begin(), totallCellList.end(), currrentCell);
      int index = std::distance(totallCellList.begin(), it1);
      for (auto cellIt = cellNeighborList[index].begin();
           cellIt != cellNeighborList[index].end();
           cellIt++)
        {
          const Particles::ParticleHandler<dim>::particle_iterator_range
particle_range = particle_handler.particles_in_cell(*cellIt); for (typename
Particles::ParticleHandler<dim>:: particle_iterator_range::iterator
partIter = particle_range.begin(); partIter != particle_range.end();
               ++partIter)
            {
              auto cPair  = std::make_pair(particleIter, partIter);
              auto cPair2 = std::make_pair(partIter, particleIter);
              auto it2 =
                std::find(contactPairs.begin(), contactPairs.end(), cPair);
              auto it3 =
                std::find(contactPairs.begin(), contactPairs.end(), cPair2);
              if (it2 == contactPairs.end())
                if (it3 == contactPairs.end())
                  if (particleIter->get_id() != partIter->get_id())
                    contactPairs.push_back(cPair);
            }
        }
    }
  return contactPairs;
}
*/

template class PPBroadSearch<2>;
template class PPBroadSearch<3>;
