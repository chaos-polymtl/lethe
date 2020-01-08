/*
 * pwcontactdetection.cpp
 *
 *  Created on: Nov 26, 2019
 *      Author: shahab
 */
#include "dem/particle_wall_contact_detection.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>

using namespace dealii;

template <int dim, int spacedim>
ParticleWallContactDetection<dim, spacedim>::ParticleWallContactDetection()
{}

template <int dim, int spacedim>
void
ParticleWallContactDetection<dim, spacedim>::boundaryCellsAndFaces(
  const Triangulation<dim, spacedim> & tr,
  std::vector<std::tuple<int,
                         typename Triangulation<dim>::active_cell_iterator,
                         int,
                         Point<dim>,
                         Point<dim>>> &boundaryCellInfo)
{
  std::vector<std::pair<int, typename Triangulation<dim>::active_cell_iterator>>
                            searchPair;
  const FE_Q<dim, spacedim> fe(1);
  QGauss<dim>               quadrature_formula(dim);
  QGauss<dim - 1>           face_quadrature_formula(dim);
  unsigned int              n_face_q_points = face_quadrature_formula.size();

  FEFaceValues<dim> fe_face_values(fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_normal_vectors);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tr.begin_active();
       cell != tr.end();
       ++cell)
    {
      for (unsigned int face_id = 0;
           face_id < GeometryInfo<dim>::faces_per_cell;
           ++face_id)
        {
          if (cell->face(face_id)->at_boundary() == true)
            {
              fe_face_values.reinit(cell, face_id);

              for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                   ++f_q_point)
                {
                  Point<dim>
                    surfNormal; // question? why doesnt it work in one line?!
                  surfNormal = fe_face_values.normal_vector(f_q_point);
                  surfNormal = -1 * surfNormal;
                  Point<dim> quadPoint = fe_face_values.quadrature_point(dim);
                  int        boundID   = cell->face(face_id)->boundary_id();
                  std::tuple cellInfo  = std::make_tuple(
                    boundID, cell, face_id, surfNormal, quadPoint);
                  std::pair cellInfoPair = std::make_pair(boundID, cell);

                  auto it = std::find(searchPair.begin(),
                                      searchPair.end(),
                                      cellInfoPair);
                  if (it == searchPair.end())
                    {
                      boundaryCellInfo.push_back(cellInfo);
                    }

                  searchPair.push_back(std::make_pair(boundID, cell));
                }
            }
        }
    }
}

template <int dim, int spacedim>
std::vector<std::tuple<
  std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
  Point<dim>,
  Point<dim>>>
ParticleWallContactDetection<dim, spacedim>::pwcontactlist(
  std::vector<std::tuple<int,
                         typename Triangulation<dim>::active_cell_iterator,
                         int,
                         Point<dim>,
                         Point<dim>>>        boundaryCellInfo,
  Particles::ParticleHandler<dim, spacedim> &particle_handler)
{
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>>>
    pwContactList;
  std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>>
    pwInfoTuple;
  // std::vector<
  //  std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>>
  //  searchPair;


  for (unsigned int i = 0; i < boundaryCellInfo.size(); ++i)
    {
      typename Triangulation<dim>::active_cell_iterator workingCell =
        std::get<1>(boundaryCellInfo[i]);
      typename Particles::ParticleHandler<dim,
                                          spacedim>::particle_iterator_range
        particle_range = particle_handler.particles_in_cell(workingCell);


      for (typename Particles::ParticleHandler<dim, spacedim>::
             particle_iterator_range::iterator partIter =
               particle_range.begin();
           partIter != particle_range.end();
           ++partIter)
        {
          std::pair pwPair =
            std::make_pair(partIter, std::get<0>(boundaryCellInfo[i]));

          //      auto it4 = std::find(searchPair.begin(), searchPair.end(),
          //      pwPair);
          //     if (it4 == searchPair.end())
          //      {
          pwInfoTuple = std::make_tuple(pwPair,
                                        std::get<3>(boundaryCellInfo[i]),
                                        std::get<4>(boundaryCellInfo[i]));
          pwContactList.push_back(pwInfoTuple);
          //     searchPair.push_back(pwPair);
          //      }
        }
    }

  return pwContactList;
}


template <int dim, int spacedim>
void
ParticleWallContactDetection<dim, spacedim>::pwFineSearch(
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>>> pwContactList,
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>,
    double,
    double,
    double,
    Point<dim>,
    double>> &   pwContactInfo,
  float          dt)
{
  std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>,
    double,
    double,
    double,
    Point<dim>,
    double>
    pwInfoTuple;
  std::vector<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>>
    pwSearchPair;
  if (!pwSearchPair.empty())
    {
      pwSearchPair.clear();
    }
  // reread this part:
  for (unsigned int i = 0; i < pwContactInfo.size(); i++)
    {
      Point<dim> partCenter =
        std::get<0>(pwContactInfo[i]).first->get_location();
      Point<dim> partVector;
      partVector = partCenter - std::get<2>(pwContactInfo[i]);
      Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactInfo[i]));

      double distance =
        ((std::get<0>(pwContactInfo[i]).first->get_properties()[2]) / 2) -
        (sqrt(projection.square()));
      if (distance > 0)
        {
          Point<dim> normVec = std::get<1>(pwContactInfo[i]);
          Point<dim> pVel    = {
            std::get<0>(pwContactInfo[i]).first->get_properties()[7],
            std::get<0>(pwContactInfo[i]).first->get_properties()[8],
            std::get<0>(pwContactInfo[i]).first->get_properties()[9]};
          Point<dim> pOmega = {
            std::get<0>(pwContactInfo[i]).first->get_properties()[16],
            std::get<0>(pwContactInfo[i]).first->get_properties()[17],
            std::get<0>(pwContactInfo[i]).first->get_properties()[18]};
          Point<dim> relVel =
            pVel +
            cross_product_3d(
              (((std::get<0>(pwContactInfo[i]).first->get_properties()[2]) /
                2) *
               pOmega),
              normVec);
          // find how to perform vec = std::get<0>(anothervec)

          double     normRelVel = relVel * normVec;
          Point<dim> relNormVel = normRelVel * normVec;
          Point<dim> relTangVel;
          relTangVel               = relVel - relNormVel;
          Point<dim> tangVec       = {0, 0, 0};
          double     relTangVelVal = relTangVel.norm();
          if (relTangVelVal != 0)
            {
              tangVec = relTangVel / relTangVelVal;
            }
          double tangRelVel = relVel * tangVec;
          double tangOverlap =
            std::get<5>(pwContactInfo[i]) + (tangRelVel * dt);

          pwInfoTuple      = std::make_tuple(std::get<0>(pwContactInfo[i]),
                                        normVec,
                                        std::get<2>(pwContactInfo[i]),
                                        distance,
                                        normRelVel,
                                        tangOverlap,
                                        tangVec,
                                        tangRelVel);
          pwContactInfo[i] = pwInfoTuple;
        }

      else
        {
          pwContactInfo.erase(pwContactInfo.begin() + i);
        }
      // find how to perform vec = std::get<0>(anothervec)
      pwSearchPair.push_back(std::get<0>(pwContactInfo[i]));
    }


  for (unsigned int i = 0; i < pwContactList.size(); i++)
    {
      Point<dim> partCenter =
        std::get<0>(pwContactList[i]).first->get_location();
      Point<dim> partVector;
      partVector = partCenter - std::get<2>(pwContactList[i]);
      Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactList[i]));
      double distance =
        ((std::get<0>(pwContactList[i]).first->get_properties()[2]) / 2) -
        (sqrt(projection.square()));

      auto it4 = std::find(pwSearchPair.begin(),
                           pwSearchPair.end(),
                           std::get<0>(pwContactList[i]));
      if (it4 == pwSearchPair.end())
        {
          if (distance > 0)
            {
              Point<dim> normVec = std::get<1>(pwContactList[i]);


              Point<dim> pVel = {
                std::get<0>(pwContactList[i]).first->get_properties()[7],
                std::get<0>(pwContactList[i]).first->get_properties()[8],
                std::get<0>(pwContactList[i]).first->get_properties()[9]};
              Point<dim> pOmega = {
                std::get<0>(pwContactList[i]).first->get_properties()[16],
                std::get<0>(pwContactList[i]).first->get_properties()[17],
                std::get<0>(pwContactList[i]).first->get_properties()[18]};
              Point<dim> relVel =
                pVel +
                cross_product_3d(
                  (((std::get<0>(pwContactList[i]).first->get_properties()[2]) /
                    2) *
                   pOmega),
                  normVec);
              // cross_prod_3d should change for 2d simulations
              double     normRelVel = relVel * normVec;
              Point<dim> relNormVel = normRelVel * normVec;
              Point<dim> relTangVel;
              relTangVel         = relVel - relNormVel;
              Point<dim> tangVec = {0, 0, 0};



              double relTangVelVal = relTangVel.norm();
              if (relTangVelVal != 0)
                {
                  tangVec = relTangVel / relTangVelVal;
                }


              double tangRelVel  = relVel * tangVec;
              double tangOverlap = 0;


              pwInfoTuple = std::make_tuple(std::get<0>(pwContactList[i]),
                                            normVec,
                                            std::get<2>(pwContactList[i]),
                                            distance,
                                            normRelVel,
                                            tangOverlap,
                                            tangVec,
                                            tangRelVel);

              pwContactInfo.push_back(pwInfoTuple);
            }
        }
    }
}


template <int dim, int spacedim>
Point<dim>
ParticleWallContactDetection<dim, spacedim>::findProjection(Point<dim> pointA,
                                                            Point<dim> pointB)
{
  Point<dim> pointC;
  pointC = ((pointA * pointB) / (pointB.square())) * pointB;

  return pointC;
}

template class ParticleWallContactDetection<3, 3>;
