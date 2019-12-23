/*
 * contactSearch.cpp
 *
 *  Created on: Oct 29, 2019
 *      Author: shahab
 */

#include "dem/contact_search.h"

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

ContactSearch::ContactSearch()
{}



std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
          std::vector<Triangulation<3>::active_cell_iterator>>
ContactSearch::findCellNeighbors(int cellNum, const Triangulation<3, 3> &tr)
{
  std::vector<std::set<Triangulation<3>::active_cell_iterator>>
                                                      cellNeighborList(cellNum);
  std::vector<Triangulation<3>::active_cell_iterator> totallCellList;

  int  iter   = 0;
  auto v_to_c = GridTools::vertex_to_cell_map(tr);
  for (Triangulation<3>::active_cell_iterator cell = tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      cellNeighborList[iter].insert(cell);
      totallCellList.push_back(cell);

      for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_cell; ++v)
        {
          for (const auto &neighbor : v_to_c[cell->vertex_index(v)])
            {
              auto it = std::find(totallCellList.begin(),
                                  totallCellList.end(),
                                  neighbor);
              if (it == totallCellList.end())
                cellNeighborList[iter].insert(neighbor);
            }
        }
      iter++;
    }
  std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
            std::vector<Triangulation<3>::active_cell_iterator>>
    contactPair = std::make_pair(cellNeighborList, totallCellList);
  return contactPair;
}



std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                      Particles::ParticleIterator<3, 3>>>
  ContactSearch::findContactPairs(
    Particles::ParticleHandler<3, 3> &                  particle_handler,
    const Triangulation<3, 3> &                         tr,
    std::vector<Triangulation<3>::active_cell_iterator> totallCellList,
    std::vector<std::set<Triangulation<3>::active_cell_iterator>>
      cellNeighborList)


// 2nd method:
{
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
      contactPairs;
  int index = 0;
  for (Triangulation<3>::active_cell_iterator cell = tr.begin_active();
       cell != tr.end();
       ++cell, ++index)
    {
      const Particles::ParticleHandler<3, 3>::particle_iterator_range
        particle_range = particle_handler.particles_in_cell(cell);


      for (auto cellIt = cellNeighborList[index].begin();
           cellIt != cellNeighborList[index].end();
           cellIt++)
        {
          const Particles::ParticleHandler<3, 3>::particle_iterator_range
            particle_range2 = particle_handler.particles_in_cell(*cellIt);

          for (typename Particles::ParticleHandler<3, 3>::
                 particle_iterator_range::iterator partIter =
                   particle_range.begin();
               partIter != particle_range.end();
               ++partIter)
            {
              for (typename Particles::ParticleHandler<3, 3>::
                     particle_iterator_range::iterator partIter2 =
                       particle_range2.begin();
                   partIter2 != particle_range2.end();
                   ++partIter2)
                {
                  auto cPair  = std::make_pair(partIter2, partIter);
                  auto cPair2 = std::make_pair(partIter, partIter2);
                  auto it2 =
                    std::find(contactPairs.begin(), contactPairs.end(), cPair);
                  auto it3 =
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


// 1st method:
/*
{
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
    contactPairs;
  if (!contactPairs.empty())
    {
      contactPairs.clear();
    }
  Triangulation<3>::active_cell_iterator currrentCell;
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
          const Particles::ParticleHandler<3, 3>::particle_iterator_range
            particle_range = particle_handler.particles_in_cell(*cellIt);
          for (typename Particles::ParticleHandler<3, 3>::
                 particle_iterator_range::iterator partIter =
                   particle_range.begin();
               partIter != particle_range.end();
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



void ContactSearch::fineSearch(
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>> contactPairs,
  dealii::Particles::ParticleHandler<3, 3> &                particle_handler,
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                   Particles::ParticleIterator<3, 3>>,
                         double,
                         Point<3>,
                         double,
                         Point<3>,
                         double,
                         double>> &                         contactInfo,
  float                                                     dt)
{
  Point<3, double> loc1, loc2;
  std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                       Particles::ParticleIterator<3, 3>>,
             double,
             Point<3>,
             double,
             Point<3>,
             double,
             double>
    infoTuple;
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
    searchPair;

  double distance;
  if (!searchPair.empty())
    {
      searchPair.clear();
    }


  for (unsigned int i = 0; i < contactInfo.size(); i++)
    {
      loc1 = (std::get<0>(contactInfo[i]).first)->get_location();
      loc2 = (std::get<0>(contactInfo[i]).second)->get_location();

      distance = ((std::get<0>(contactInfo[i]).first->get_properties()[2] +
                   std::get<0>(contactInfo[i]).second->get_properties()[2]) /
                  2) -
                 loc1.distance(loc2);
      if (distance > 0)
        {
          Point<3> contactVector;
          contactVector    = (loc2 - loc1);
          Point<3> normVec = contactVector / contactVector.norm();

          Point<3> part1Vel = {
            std::get<0>(contactInfo[i]).first->get_properties()[7],
            std::get<0>(contactInfo[i]).first->get_properties()[8],
            std::get<0>(contactInfo[i]).first->get_properties()[9]};
          Point<3> part2Vel = {
            std::get<0>(contactInfo[i]).second->get_properties()[7],
            std::get<0>(contactInfo[i]).second->get_properties()[8],
            std::get<0>(contactInfo[i]).second->get_properties()[9]};

          Point<3> part1AngVel = {
            std::get<0>(contactInfo[i]).first->get_properties()[16],
            std::get<0>(contactInfo[i]).first->get_properties()[17],
            std::get<0>(contactInfo[i]).first->get_properties()[18]};

          Point<3> part2AngVel = {
            std::get<0>(contactInfo[i]).second->get_properties()[16],
            std::get<0>(contactInfo[i]).second->get_properties()[17],
            std::get<0>(contactInfo[i]).second->get_properties()[18]};

          // ************************************************
          // in tempelate <3> <2> take care
          Point<3> relVel;
          relVel =
            (part1Vel - part2Vel) +
            (cross_product_3d(
              (((std::get<0>(contactInfo[i]).first->get_properties()[2] / 2.0) *
                part1AngVel) +
               ((std::get<0>(contactInfo[i]).second->get_properties()[2] /
                 2.0) *
                part2AngVel)),
              normVec));

          double   normRelVel = relVel * normVec;
          Point<3> relNormVel = normRelVel * normVec;
          Point<3> relTangVel;
          relTangVel       = relVel - relNormVel;
          Point<3> tangVec = {0, 0, 0};

          double relTangVelVal = relTangVel.norm();
          if (relTangVelVal != 0)
            {
              tangVec = relTangVel / relTangVelVal;
            }

          double tangRelVel  = relVel * tangVec;
          double tangOverlap = std::get<6>(contactInfo[i]) + (tangRelVel * dt);

          infoTuple      = std::make_tuple(std::get<0>(contactInfo[i]),
                                      distance,
                                      normVec,
                                      normRelVel,
                                      tangVec,
                                      tangRelVel,
                                      tangOverlap);
          contactInfo[i] = infoTuple;
        }
      else
        {
          contactInfo.erase(contactInfo.begin() + i);
        }
      searchPair.push_back(std::get<0>(contactInfo[i]));
    }



  for (unsigned int i = 0; i < contactPairs.size(); i++)
    {
      loc1 = (contactPairs[i].first)->get_location();
      loc2 = (contactPairs[i].second)->get_location();

      distance = ((contactPairs[i].first->get_properties()[2] +
                   contactPairs[i].second->get_properties()[2]) /
                  2) -
                 loc1.distance(loc2);

      auto it4 =
        std::find(searchPair.begin(), searchPair.end(), contactPairs[i]);
      if (it4 == searchPair.end())
        {
          if (distance > 0)
            {
              Point<3> contactVector;
              contactVector    = (loc2 - loc1);
              Point<3> normVec = contactVector / contactVector.norm();

              Point<3> part1Vel = {contactPairs[i].first->get_properties()[7],
                                   contactPairs[i].first->get_properties()[8],
                                   contactPairs[i].first->get_properties()[9]};
              Point<3> part2Vel = {contactPairs[i].second->get_properties()[7],
                                   contactPairs[i].second->get_properties()[8],
                                   contactPairs[i].second->get_properties()[9]};

              Point<3> part1AngVel = {
                contactPairs[i].first->get_properties()[16],
                contactPairs[i].first->get_properties()[17],
                contactPairs[i].first->get_properties()[18]};

              Point<3> part2AngVel = {
                contactPairs[i].second->get_properties()[16],
                contactPairs[i].second->get_properties()[17],
                contactPairs[i].second->get_properties()[18]};

              // ************************************************
              // in tempelate <3> <2> take care
              Point<3> relVel;
              relVel = (part1Vel - part2Vel) +
                       (cross_product_3d(
                         (((contactPairs[i].first->get_properties()[2] / 2.0) *
                           part1AngVel) +
                          ((contactPairs[i].second->get_properties()[2] / 2.0) *
                           part2AngVel)),
                         normVec));



              double   normRelVel = relVel * normVec;
              Point<3> relNormVel = normRelVel * normVec;
              Point<3> relTangVel;
              relTangVel = relVel - relNormVel;

              Point<3> tangVec       = {0, 0, 0};
              double   relTangVelVal = relTangVel.norm();
              if (relTangVelVal != 0)
                {
                  tangVec = relTangVel / relTangVelVal;
                }
              double tangRelVel  = relVel * tangVec;
              double tangOverlap = 0;

              infoTuple =
                std::make_tuple(std::make_pair(contactPairs[i].first,
                                               contactPairs[i].second),
                                distance,
                                normVec,
                                normRelVel,
                                tangVec,
                                tangRelVel,
                                tangOverlap);
              contactInfo.push_back(infoTuple);
            }
        }
    }
}
