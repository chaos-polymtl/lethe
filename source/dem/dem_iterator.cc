/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "dem/dem_iterator.h"

#include <deal.II/particles/particle_handler.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include "dem/contact_force.h"
#include "dem/contact_search.h"
#include "dem/integration.h"
#include "dem/parameters_dem.h"
#include "dem/particle_wall_contact_detection.h"
#include "dem/particle_wall_contact_force.h"
#include "dem/visualization.h"
#include "dem/write_vtu.h"
#include <deal.II/particles/property_pool.h>

#include <math.h>


using namespace dealii;

DEM_iterator::DEM_iterator()
{}

void
  DEM_iterator::forceReinit(Particles::ParticleHandler<3, 3> &particle_handler)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      particle->get_properties()[13] = 0;
      particle->get_properties()[14] = 0;
      particle->get_properties()[15] = 0;
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

void
DEM_iterator::engine(
  int &                             nPart,
  Particles::ParticleHandler<3, 3> &particle_handler,
  const Triangulation<3, 3> &       tr,
  int &                             step,
  float &                           time,
  ParametersDEM<3>                  DEMparam,
  std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
            std::vector<Triangulation<3>::active_cell_iterator>> cellNeighbor,
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                   Particles::ParticleIterator<3, 3>>,
                         double,
                         Point<3>,
                         double,
                         Point<3>,
                         double,
                         double>> &                              contactInfo,
  std::vector<std::tuple<int,
                         Triangulation<3>::active_cell_iterator,
                         int,
                         Point<3>,
                         Point<3>>> boundaryCellInfo,
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                         Point<3>,
                         Point<3>,
                         double,
                         double,
                         double,
                         Point<3>,
                         double>> & pwContactInfo, std::vector<std::tuple<std::string, int>> properties, Particles::PropertyPool &propPool)
{
  // moving walls


  // check simulation boundaries
  // checkSimBound(particle_handler, readInput);

    //insertion
    if (fmod(step, DEMparam.insertionInfo.insertFrequncy) == 1)
      {
    if (step < DEMparam.insertionInfo.tInsertion)
      {
        if (nPart < DEMparam.simulationControl.nTotal) // number < total number
          {
  ParticleInsertion ins1(DEMparam);
                ins1.uniformInsertion(
                  particle_handler, tr, DEMparam, nPart, propPool);
              }
          }
    }


  // contact search
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
    contactPairs;

  // force reinitilization
  forceReinit(particle_handler);

  ContactSearch cs;
  // if (fmod(step,10) == 1)
  //	{
  contactPairs = cs.findContactPairs(particle_handler,
                                     tr,
                                     cellNeighbor.second,
                                     cellNeighbor.first);
  //	}

  cs.fineSearch(contactPairs,
                particle_handler,
                contactInfo,
                DEMparam.simulationControl.dt);

  // contact force
  ContactForce cf;
  cf.linearCF(contactInfo, particle_handler, DEMparam);

  // p-w contact detection:
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                         Point<3>,
                         Point<3>>>
    pwContactList;


  ParticleWallContactDetection pw;
  // if (fmod(step,10) == 1)
  //{
  pwContactList = pw.pwcontactlist(boundaryCellInfo, particle_handler);
  //}


  pw.pwFineSearch(pwContactList,
                  particle_handler,
                  pwContactInfo,
                  DEMparam.simulationControl.dt);


  // p-w contact force:
  ParticleWallContactForce pwcf;
  pwcf.pwLinearCF(pwContactInfo, DEMparam);


  // Integration
  Integration Integ1;
  Integ1.eulerIntegration(particle_handler, DEMparam);
  //Integ1.rk2Integration(particle_handler, DEMparam);



  //visualization
  if (fmod(step, DEMparam.simulationControl.writeFrequency) == 1)
    {
      Visualization visObj;
      visObj.build_patches(particle_handler,
                           DEMparam.outputProperties.numFields,
                           DEMparam.outputProperties.numProperties,
                           properties);
      WriteVTU writObj;
      writObj.write_master_files(visObj);
      writObj.writeVTUFiles(visObj, step, time);
    }

  //print iteration
  if (fmod(step,1000) == 0)
  {
            std::cout << "Step "<< step << std::endl;
  }

  // update:
  particle_handler.sort_particles_into_subdomains_and_cells();
  step = step + 1;
  time = step * DEMparam.simulationControl.dt;
}
