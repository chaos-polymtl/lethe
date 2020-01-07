/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "dem/dem_iterator.h"

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <math.h>

#include "dem/parameters_dem.h"
#include "dem/visualization.h"
#include "dem/write_vtu.h"


using namespace dealii;

template <int dim, int spacedim>
DEM_iterator<dim, spacedim>::DEM_iterator()
{}

template <int dim, int spacedim>
void
DEM_iterator<dim, spacedim>::forceReinit(
  Particles::ParticleHandler<dim, spacedim> &particle_handler)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      particle->get_properties()[13] = 0;
      particle->get_properties()[14] = 0;
      particle->get_properties()[15] = 0;
      /*
      particle->get_properties()[21] = 0;
      particle->get_properties()[22] = 0;
      particle->get_properties()[23] = 0;
      */

    }
}

/*
void DEM_iterator::checkSimBound(
  Particles::ParticleHandler<3, 3> &particle_handler,
  ReadInputScript                   readInput)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      if (particle->get_properties()[4] < readInput.x_min ||
          particle->get_properties()[4] > readInput.x_max ||
          particle->get_properties()[5] < readInput.y_min ||
          particle->get_properties()[5] > readInput.y_max ||
          particle->get_properties()[6] < readInput.z_min ||
          particle->get_properties()[6] > readInput.z_max)
        {
          particle_handler.remove_particle(particle);
        }
    }
}
*/


template <int dim, int spacedim>
void
DEM_iterator<dim, spacedim>::engine(
  int &                                      nPart,
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const Triangulation<dim, spacedim> &       tr,
  int &                                      step,
  float &                                    time,
  std::pair<
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>,
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
                                      cellNeighbor,
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>,
                                   Particles::ParticleIterator<dim, spacedim>>,
                         double,
                         Point<dim>,
                         double,
                         Point<dim>,
                         double,
                         double>> &   contactInfo,
  std::vector<std::tuple<int,
                         typename Triangulation<dim>::active_cell_iterator,
                         int,
                         Point<dim>,
                         Point<dim>>> boundaryCellInfo,
  std::vector<
    std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>, int>,
               Point<dim>,
               Point<dim>,
               double,
               double,
               double,
               Point<dim>,
               double>> &                     pwContactInfo,
  std::vector<std::tuple<std::string, int>>   properties,
        Particles::PropertyPool &                   propPool,
  ContactSearch<dim, spacedim>                cs,
  ParticleWallContactDetection<dim, spacedim> pw,
  ContactForce<dim, spacedim>                 cf,
  ParticleWallContactForce<dim, spacedim>     pwcf,
  Integration<dim, spacedim>                  Integ1, int numberOfSteps, double dt, int nTotal, int writeFreq, Point<dim> g, double dp, int rhop, int Yp, int Yw, float vp, float vw, float ep, float ew, float mup, float muw, float murp, float murw, int tInsertion, int nInsert, int insertFrequency, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max, int numFields, int numProperties)
{
  // moving walls


  // insertion
  if (fmod(step, insertFrequency) == 1)
    {
      if (step < tInsertion)
        {
          if (nPart <
              nTotal) // number < total number
            {
              ParticleInsertion<dim, spacedim> ins1(x_min, y_min,z_min, x_max, y_max, z_max, dp, nInsert);
              ins1.nonUniformInsertion(
                particle_handler, tr, nPart, propPool, x_min, y_min, z_min, x_max, y_max, z_max, dp, nInsert,rhop, g);
            }
        }
    }


  // contact search
  std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                        Particles::ParticleIterator<dim, spacedim>>>
    contactPairs;

  // force reinitilization
  forceReinit(particle_handler);

  if (fmod(step,10) == 1)
        {
  contactPairs = cs.findContactPairs(particle_handler,
                                     tr,
                                     cellNeighbor.first);
    }

  cs.fineSearch(contactPairs, contactInfo, dt);

  // contact force
  cf.nonLinearCF(contactInfo, Yp, vp, mup, murp);

  // p-w contact detection:
  std::vector<
    std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>, int>,
               Point<dim>,
               Point<dim>>>
    pwContactList;


 if (fmod(step,10) == 1)
 {
  pwContactList = pw.pwcontactlist(boundaryCellInfo, particle_handler);
 }


  pw.pwFineSearch(pwContactList, pwContactInfo, dt);


  // p-w contact force:
  pwcf.pwNonLinearCF(pwContactInfo, vp, Yp, vw, Yw, muw, murw);


  // Integration
   Integ1.eulerIntegration(particle_handler, g, dt);
  //Integ1.rk2Integration(particle_handler, g, dt);


  // visualization
  if (fmod(step, writeFreq) == 1)
    {
      Visualization<dim, spacedim> visObj;
      visObj.build_patches(particle_handler,
                           numFields,
                           numProperties,
                           properties);
      WriteVTU<dim, spacedim> writObj;
      writObj.write_master_files(visObj);
      writObj.writeVTUFiles(visObj, step, time);
    }

  // print iteration
  if (fmod(step, 1000) == 0)
    {
      std::cout << "Step " << step << std::endl;
    }

  // update:
  particle_handler.sort_particles_into_subdomains_and_cells();
  step = step + 1;
  time = step * dt;
}

template class DEM_iterator<3, 3>;
