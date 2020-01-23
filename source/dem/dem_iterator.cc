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

#include <chrono>

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
  Integration<dim, spacedim>                  Integ1,
  int                                         numberOfSteps,
  double                                      dt,
  int                                         nTotal,
  int                                         writeFreq,
  Point<dim>                                  g,
  double                                      dp,
  int                                         rhop,
  int                                         Yp,
  int                                         Yw,
  float                                       vp,
  float                                       vw,
  float                                       ep,
  float                                       ew,
  float                                       mup,
  float                                       muw,
  float                                       murp,
  float                                       murw,
  int                                         tInsertion,
  int                                         nInsert,
  int                                         insertFrequency,
  float                                       x_min,
  float                                       y_min,
  float                                       z_min,
  float                                       x_max,
  float                                       y_max,
  float                                       z_max,
  int                                         numFields,
  int                                         numProperties)
{
  // moving walls

  auto t1 = std::chrono::high_resolution_clock::now();
  // insertion
  if (fmod(step, insertFrequency) == 1)
    {
      if (step < tInsertion)
        {
          if (nPart < nTotal) // number < total number
            {
              ParticleInsertion<dim, spacedim> ins1(
                x_min, y_min, z_min, x_max, y_max, z_max, dp, nInsert);
              ins1.nonUniformInsertion(particle_handler,
                                       tr,
                                       nPart,
                                       propPool,
                                       x_min,
                                       y_min,
                                       z_min,
                                       x_max,
                                       y_max,
                                       z_max,
                                       dp,
                                       nInsert,
                                       rhop,
                                       g);
            }
        }
    }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration_Insertion =
    std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  // contact search
  std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                        Particles::ParticleIterator<dim, spacedim>>>
    contactPairs;

  // force reinitilization
  forceReinit(particle_handler);

  auto t3 = std::chrono::high_resolution_clock::now();
  // if (fmod(step,10) == 1)
  //      {
  contactPairs = cs.findContactPairs(particle_handler, tr, cellNeighbor.first);
  //  }
  auto t4 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactPairs =
    std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();

  auto t5 = std::chrono::high_resolution_clock::now();
  cs.fineSearch(contactPairs, contactInfo, dt);
  auto t6 = std::chrono::high_resolution_clock::now();
  auto duration_PPFineSearch =
    std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();

  auto t7 = std::chrono::high_resolution_clock::now();
  // contact force
  cf.nonLinearCF(contactInfo, Yp, vp, mup, murp);
  // cf.linearCF(contactInfo, Yp, vp, ep, mup, murp);
  auto t8 = std::chrono::high_resolution_clock::now();
  auto duration_PPContactForce =
    std::chrono::duration_cast<std::chrono::microseconds>(t8 - t7).count();

  // p-w contact detection:
  std::vector<
    std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>, int>,
               Point<dim>,
               Point<dim>>>
    pwContactList;


  auto t9 = std::chrono::high_resolution_clock::now();
  // if (fmod(step,10) == 1)
  //{
  pwContactList = pw.pwcontactlist(boundaryCellInfo, particle_handler);
  // }
  auto t10 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactPairs =
    std::chrono::duration_cast<std::chrono::microseconds>(t10 - t9).count();

  auto t11 = std::chrono::high_resolution_clock::now();
  pw.pwFineSearch(pwContactList, pwContactInfo, dt);
  auto t12 = std::chrono::high_resolution_clock::now();
  auto duration_PWFineSearch =
    std::chrono::duration_cast<std::chrono::microseconds>(t12 - t11).count();

  auto t13 = std::chrono::high_resolution_clock::now();
  // p-w contact force:
  pwcf.pwNonLinearCF(pwContactInfo, vp, Yp, vw, Yw, muw, murw);
  // pwcf.pwLinearCF(pwContactInfo, vp, Yp, vw, Yw, ew, muw, murw);
  auto t14 = std::chrono::high_resolution_clock::now();
  auto duration_PWContactForce =
    std::chrono::duration_cast<std::chrono::microseconds>(t14 - t13).count();

  auto t15 = std::chrono::high_resolution_clock::now();
  // Integration
  // Integ1.eulerIntegration(particle_handler, g, dt);
  // Integ1.rk2Integration(particle_handler, g, dt);
  Integ1.velocity_verlet_integration(particle_handler, g, dt);
  auto t16 = std::chrono::high_resolution_clock::now();
  auto duration_Integration =
    std::chrono::duration_cast<std::chrono::microseconds>(t16 - t15).count();

  auto t17 = std::chrono::high_resolution_clock::now();
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
  auto t18 = std::chrono::high_resolution_clock::now();
  auto duration_Visualization =
    std::chrono::duration_cast<std::chrono::microseconds>(t18 - t17).count();

  // print iteration
  if (fmod(step, 1000) == 1)
    {
      std::cout << "Step " << step << std::endl;
      std::cout << "CPU time of insertion is: " << duration_Insertion
                << " micros" << std::endl;
      std::cout << "CPU time of P-P borad search is: "
                << duration_PPContactPairs << " micros" << std::endl;
      std::cout << "CPU time of P-P fine search is: " << duration_PPFineSearch
                << " micros" << std::endl;
      std::cout << "CPU time of P-P contact force is: "
                << duration_PPContactForce << " micros" << std::endl;
      std::cout << "CPU time of P-W borad search is: "
                << duration_PWContactPairs << " micros" << std::endl;
      std::cout << "CPU time of P-W fine search is: " << duration_PWFineSearch
                << " micros" << std::endl;
      std::cout << "CPU time of P-W contact force is: "
                << duration_PWContactForce << " micros" << std::endl;
      std::cout << "CPU time of integration is: " << duration_Integration
                << " micros" << std::endl;
      std::cout << "CPU time of visualization is: " << duration_Visualization
                << " micros" << std::endl;
      std::cout
        << "-----------------------------------------------------------------------"
        << std::endl;
    }

  // update:
  particle_handler.sort_particles_into_subdomains_and_cells();
  step = step + 1;
  time = step * dt;
}

template class DEM_iterator<3, 3>;
